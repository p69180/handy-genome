import logging

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))
workflow = importlib.import_module('.'.join([top_package_name, 'workflow']))


DEFAULT_NUM_SUBMITTED = 400


def argument_parser(cmdargs):
	def sanity_check(args):
		pass

	parser_dict = workflow.init_parser()
	workflow.add_indir_arg(
			parser_dict['required'], 
			required = True,
			help = f'Directory where slurm job scripts exist.',
			)
	workflow.add_checkint_arg(parser_dict['optional'], required = False) # --slurm-check-interval, args.intv_check
	workflow.add_submitint_arg(parser_dict['optional'], required = False) # --slurm-submit-interval, args.intv_submit
	workflow.add_logging_args(parser_dict)

	parser_dict['optional'].add_argument(
			'-n', required = False, 
			default = DEFAULT_NUM_SUBMITTED, type = int, 
			help = f'The number of jobs to submitted simultaneously. Default: {DEFAULT_NUM_SUBMITTED}',
			)
	parser_dict['flag'].add_argument(
			'--nowait', action = 'store_true',
			help = f'If set, this program will not wait for the submitted jobs to finish.',
			)

	args = parser_dict['main'].parse_args(cmdargs)
	sanity_check(args)

	return args


def main(cmdargs):
	args = argument_parser(cmdargs)
	logger = workflow.get_logger(
			name = 'sbatch-helper', 
			stderr = (not args.silent),
			filename = args.log,
			append = False,
			)

	jobscript_path_list = common.listdir(args.indir_path)
	joblist = workflow.JobList(jobscript_path_list, logger = logger)
	joblist.submit(intv_submit = args.intv_submit)
	if not args.nowait:
		joblist.wait(intv_check = args.intv_check, edit_log_suffix = True)

	logger.info('ALL FINISHED SUCCESSFULLY')
