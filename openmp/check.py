#!/usr/bin/python

if __name__ == '__main__':
	with open("../OutputData/RlxMthd_v1.0_300.answere","r") as answere_file:
		answere_lines = answere_file.read()
		with open('../OutputData/RlxMthd_v1.0_300.dat','r') as output_file:
			output_lines = output_file.read()
			print "Test Success !" if output_lines == answere_lines else "Test Fail !"