import os

work_dir = "/home/users/niv/chimeric_help/results/U_motifs"
files = os.listdir(work_dir)

for name in files:

	file_path = os.path.join(work_dir, name)

	print file_path
	if file_path.endswith(".txt"):

		os.system("python generate_logo.py -f %s" % file_path)