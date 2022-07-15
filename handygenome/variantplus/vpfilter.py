import re
import functools

import pysam
import numpy as np

import importlib

top_package_name = __name__.split(".")[0]
common = importlib.import_module(".".join([top_package_name, "common"]))


FORMAT_FILTER_META = {
    "ID": "filter",
    "Type": "String",
    "Number": ".",
    "Description": (f"Works as FILTER column. The variant is true one if only 'PASS' exists."),
}


def add_format_filter_meta(vcfheader):
    vcfheader.add_meta(key="FORMAT", items=FORMAT_FILTER_META.items())


class VariantPlusFilter:
    def __repr__(self):
        return f"<{self.__class__.__name__} ({self.params})>"


class SamplewiseFilter(VariantPlusFilter):
    @classmethod
    def add_filter(cls, vp, sampleid):
        vp.add_sample_filter(sampleid, cls.__name__)

    def apply(self, vp, sampleid):
        mask = self.check(vp, sampleid)
        if not mask:
            self.add_filter(vp, sampleid)

    def check_samples(self, vp, sampleids=None, alleleindex=1):
        if sampleids is None:
            sampleids = sorted(vp.vr.samples.keys())
        
        return np.array([self.check(vp=vp, sampleid=sampleid, alleleindex=alleleindex)
                         for sampleid in sampleids])

    def show_result(self, vp, sampleid, alleleindex=1):
        mask = self.check(vp, sampleid, alleleindex)
        print(mask, self, sep=(' ' * 4))

    def _check_meanvalue_difference(
        self, readstats_key, vp, sampleid, alleleindex, cutoff
    ):
        readstats = vp.readstats_dict[sampleid]
        other_alleleindexes = vp.get_other_alleleindexes(alleleindex)

        target_value = readstats[readstats_key][alleleindex]
        other_value = readstats.get_alleleindexes_mean(
            readstats_key, other_alleleindexes
        )

        if np.isnan([target_value, other_value]).any():
            return True
        else:
            return (target_value - other_value) >= cutoff

    def _check_absvalue(self, readstats_key, vp, sampleid, alleleindex, cutoff):
        readstats = vp.readstats_dict[sampleid]
        target_value = readstats[readstats_key][alleleindex]

        if np.isnan(target_value):
            return True
        else:
            return target_value >= cutoff


class NonSamplewiseFilter(VariantPlusFilter):
    def show_result(self, vp):
        mask = self.check(vp)
        print(mask, self, sep=(' ' * 4))


#######################################


class PopfreqFilter(NonSamplewiseFilter):
    def __init__(
        self, popnames=("GnomAD", "1000Genomes", "KOREAN", "Korea1K"), cutoff=0.01
    ):
        self.params = {
            "popnames": popnames,
            "cutoff": cutoff,
        }

    def check(self, vp):
        return all(
            vp.get_popfreq(popname) <= self.params["cutoff"]
            for popname in self.params["popnames"]
        )


class PonFilter(SamplewiseFilter):
    @common.get_deco_arg_choices({"mode": ("mean", "max", "median")})
    def __init__(
        self,
        readstats_dict,
        deviation_cutoff=5,
        subset_num_cutoff=3,
        subset_fraction=0.5,
        germline_vaf_cutoff=0.25,
        germline_sample_ratio_cutoff=0.4,
        mode="mean",
    ):
        self.readstats_dict = readstats_dict
        self.params = {
            "deviation_cutoff": deviation_cutoff,
            "subset_num_cutoff": subset_num_cutoff,
            "subset_fraction": subset_fraction,
            "germline_vaf_cutoff": germline_vaf_cutoff,
            "germline_sample_ratio_cutoff": germline_sample_ratio_cutoff,
            "mode": mode,
        }

        self.sampleids = tuple(sorted(self.readstats_dict.keys()))
        self.nsample = len(self.readstats_dict)

    def __repr__(self):
        infostr = ", ".join(
            [
                f"sample number: {self.nsample}",
                f"params: {self.params}",
                f"sample IDs: {self.sampleids}",
            ]
        )
        return f"<PonFilter ({infostr})>"

    @functools.cache
    def get_vaf_dict(self, alleleindex):
        vaf_dict = dict()
        for sampleid in self.sampleids:
            readstats = self.readstats_dict[sampleid]
            denom = readstats.get_total_rppcount()
            if denom == 0:
                vaf = np.nan
            else:
                numer = readstats["rppcounts"][alleleindex]
                vaf = numer / denom
            vaf_dict[sampleid] = vaf

        return vaf_dict

    @functools.cache
    def get_vafs(self, alleleindex):
        vaf_dict = self.get_vaf_dict(alleleindex)
        vafs = [vaf_dict[key] for key in self.sampleids]

        return np.array(vafs)

    @functools.cache
    def get_noise_summary(self, alleleindex):
        altcnts = np.array(
            [
                self.readstats_dict[sampleid]["rppcounts"][alleleindex]
                for sampleid in self.sampleids
            ]
        )
        altcnts_nonzero = altcnts[altcnts > 0]
        subset_num = int(len(altcnts_nonzero) * self.params["subset_fraction"])

        if subset_num >= self.params["subset_num_cutoff"]:
            subset = sorted(altcnts_nonzero)[:subset_num]
            noise_summary = getattr(np, self.params["mode"])(subset)
        else:
            noise_summary = None

        return noise_summary

    def check_greater_than_noise(self, vp, sampleid, alleleindex):
        allele_count = vp.readstats_dict[sampleid]["rppcounts"][alleleindex]
        noise_summary = self.get_noise_summary(alleleindex)
        if noise_summary is None:
            return True
        else:
            deviation = allele_count / noise_summary
            return deviation >= self.params["deviation_cutoff"]

    @functools.cache
    def get_noise_summary_vaf(self, alleleindex):
        vafs = self.get_vafs(alleleindex)
        vafs_nonzero = vafs[vafs > 0]
        subset_num = int(len(vafs_nonzero) * self.params["subset_fraction"])

        if subset_num >= self.params["subset_num_cutoff"]:
            subset = sorted(vafs_nonzero)[:subset_num]
            noise_summary = getattr(np, self.params["mode"])(subset)
        else:
            noise_summary = None

        return noise_summary

    def check_greater_than_noise_vaf(self, vp, sampleid, alleleindex):
        vaf = self.get_vaf_dict(alleleindex)[sampleid]
        noise_summary = self.get_noise_summary_vaf(alleleindex)
        if noise_summary is None:
            return True
        else:
            deviation = vaf / noise_summary
            return deviation >= self.params["deviation_cutoff"]

    @functools.cache
    def check_germline_pattern(self, alleleindex):
        vafs = self.get_vafs(alleleindex)
        n_highvaf = (vafs >= self.params["germline_vaf_cutoff"]).sum()
        highvaf_ratio = n_highvaf / self.nsample

        return highvaf_ratio >= self.params["germline_sample_ratio_cutoff"]

    @functools.cache
    def check_germline_pattern_from_self(self, alleleindex):
        vafs = self.get_vafs(alleleindex)
        noise_summary = self.get_noise_summary_vaf(alleleindex)
        if noise_summary is None:
            return False
        else:
            n_highvaf = (vafs >= noise_summary).sum()
            highvaf_ratio = n_highvaf / self.nsample
            return highvaf_ratio >= self.params["germline_sample_ratio_cutoff"]

    def check(self, vp, sampleid, alleleindex=1):
        #is_germline_pattern = self.check_germline_pattern(alleleindex)
        is_germline_pattern = self.check_germline_pattern_from_self(alleleindex)
        #is_gt_noise = self.check_greater_than_noise(vp, sampleid, alleleindex)
        is_gt_noise = self.check_greater_than_noise_vaf(vp, sampleid, alleleindex)

        return ((not is_germline_pattern) and is_gt_noise)


class DiffMeanBQFilter(SamplewiseFilter):
    def __init__(self, cutoff=-5):
        self.params = {
            "cutoff": cutoff,
        }

    def check(self, vp, sampleid, alleleindex=1):
        return self._check_meanvalue_difference(
            readstats_key="mean_BQs",
            vp=vp,
            sampleid=sampleid,
            alleleindex=alleleindex,
            cutoff=self.params["cutoff"],
        )


class AbsMeanBQFilter(SamplewiseFilter):
    def __init__(self, cutoff=20):
        self.params = {
            "cutoff": cutoff,
        }

    def check(self, vp, sampleid, alleleindex=1):
        return self._check_absvalue(
            readstats_key="mean_BQs",
            vp=vp,
            sampleid=sampleid,
            alleleindex=alleleindex,
            cutoff=self.params["cutoff"],
        )


class DiffMeanMQFilter(SamplewiseFilter):
    def __init__(self, cutoff=-15):
        self.params = {
            "cutoff": cutoff,
        }

    def check(self, vp, sampleid, alleleindex=1):
        return self._check_meanvalue_difference(
            readstats_key="mean_MQs",
            vp=vp,
            sampleid=sampleid,
            alleleindex=alleleindex,
            cutoff=self.params["cutoff"],
        )


class AbsMeanMQFilter(SamplewiseFilter):
    def __init__(self, cutoff=40):
        self.params = {
            "cutoff": cutoff,
        }

    def check(self, vp, sampleid, alleleindex=1):
        return self._check_absvalue(
            readstats_key="mean_MQs",
            vp=vp,
            sampleid=sampleid,
            alleleindex=alleleindex,
            cutoff=self.params["cutoff"],
        )


class ClipoverlapFilter(SamplewiseFilter):
    def __init__(self, cutoff=1):
        self.params = {
            "cutoff": cutoff,
        }

    def check(self, vp, sampleid, alleleindex=1):
        readstats = vp.readstats_dict[sampleid]
        clipoverlap_count = readstats["rppcounts"]["softclip_overlap"]
        target_count = readstats["rppcounts"][alleleindex]

        if clipoverlap_count == 0:
            return True
        else:
            ratio = target_count / clipoverlap_count
            return ratio > self.params["cutoff"]


class VarposUniformFilter(SamplewiseFilter):
    def __init__(self, cutoff=0.05):
        self.params = {
            "cutoff": cutoff,
        }

    def check(self, vp, sampleid, alleleindex=1):
        readstats = vp.readstats_dict[sampleid]
        pval = readstats["varpos_uniform_pvalues"][alleleindex]
        if np.isnan(pval):
            return True
        else:
            return pval >= self.params["cutoff"]


class VarposValueFilter(SamplewiseFilter):
    def __init__(self, cutoff=0.3):
        self.params = {
            "cutoff": cutoff,
        }

    def check(self, vp, sampleid, alleleindex=1):
        return self._check_absvalue(
            readstats_key="mean_varpos_fractions",
            vp=vp,
            sampleid=sampleid,
            alleleindex=alleleindex,
            cutoff=self.params["cutoff"],
        )


class ReadcountFilter(SamplewiseFilter):
    def __init__(self, cutoff=2):
        self.params = {
            "cutoff": cutoff,
        }

    def check(self, vp, sampleid, alleleindex=1):
        return self._check_absvalue(
            readstats_key="rppcounts",
            vp=vp,
            sampleid=sampleid,
            alleleindex=alleleindex,
            cutoff=self.params["cutoff"],
        )


class OthercountRatioFilter(SamplewiseFilter):
    def __init__(self, cutoff=1.5):
        self.params = {
            "cutoff": cutoff,
        }

    def check(self, vp, sampleid, alleleindex=1):
        readstats = vp.readstats_dict[sampleid]
        target_count = readstats["rppcounts"][alleleindex]
        other_count = readstats["rppcounts"][-1]

        if other_count == 0:
            return True
        else:
            ratio = target_count / other_count
            return ratio > self.params["cutoff"]


############################


def get_transcript_subtype_filter(key):
    def vpfilter(vp):
        return any(
            feature["transcript_subtype_flags"][key]
            for feature in vp.annotdb.transcript_canon_ovlp.values()
        )

    return vpfilter


CODING_GENE_INVOLVED = get_transcript_subtype_filter("is_coding")


def get_consequence_filter(key):
    def vpfilter(vp):
        return any(
            feature["consequence_flags"][key]
            for feature in vp.annotdb.transcript_canon_ovlp.values()
        )

    return vpfilter


PROTEIN_ALTERING = get_consequence_filter("is_protein_altering")
NOT_PROTEIN_ALTERING = get_consequence_filter("is_not_protein_altering")


def get_genename_filter(gene_list, canonical_only=True):
    gene_set = set(gene_list)
    if canonical_only:
        def vpfilter(vp):
            vp_genes = set(vp.get_gene_names(canonical_only=True))
            return not vp_genes.isdisjoint(gene_set)
    else:
        def vpfilter(vp):
            vp_genes = set(vp.get_gene_names(canonical_only=False))
            return not vp_genes.isdisjoint(gene_set)

    return vpfilter
