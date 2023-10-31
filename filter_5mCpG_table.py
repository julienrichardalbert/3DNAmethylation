# python ./filter_build_on_me.py > input_file_filtered.txt
# 0 = 1
# x = 24
import sys
import numpy

firstarg=sys.argv[1]

# number of first data column (usually COUNT) minus one
first_data_column = 13
# number of total columns minus two
last_data_column = 20
column_indices = numpy.arange(first_data_column, last_data_column, 3)

# defaults
cg_count_filter = 1
range_filter    = 1

# import file
with open(firstarg, 'r') as txtfile:
# pop off first line from file and store as variable "header"
	header = txtfile.readline()
	# only keep first "first_data_column" columns and every 4 columns after that
	newHeader = header.split("	")[0:first_data_column]
	for index in column_indices:
		newHeader.append(header.split("	")[index+1])
	# print the list as a tab delimited row of values
	print('\t'.join(map(str,newHeader)))
#	print newHeader

	for row in txtfile:
		listToPrint=[]
		# keep first "first_data_column" which we will not filter on
		for first_cols in range(0,first_data_column):
			listToPrint.append(row.split("	")[first_cols])

		# set up variables to filter on
		for index in column_indices:
			cg_count       = int(row.split("	")[index])
			mean_meth      =     row.split("	")[index+1]
			coverage_range = int(row.split("	")[index+2])

			# get dataset name, then set filter values accordingly
			dataset = header.split("	")[index]
#			if "x5" in dataset:
#				cpg_count_filter = 2
#			if "x1" in dataset:
#				cpg_count_filter = 5
#			if "Wang2014" in dataset:
#				coverage_range = 500
#			if "Kobayashi2013" in dataset:
#				coverage_range = 500
#			if "Wang2019" in dataset:
#				coverage_range = 500

			# actual filtering happens here
			if cg_count >= cg_count_filter and coverage_range >= range_filter:
				newVal = mean_meth
			else:
				newVal = "NaN"
			listToPrint.append(newVal)

		# print the list as a tab delimited row of values
		print('\t'.join(map(str,listToPrint)))
#		print listToPrint
