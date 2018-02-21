import re

pattern = ' -308 {3,7}-309 {3,7}151 '

with open("log.txt", "r") as ins:
	line_no = 0

	for line in ins:
		line_no += 1
		if re.search(pattern, line):
			print('line_no=', line_no)
			print(line)