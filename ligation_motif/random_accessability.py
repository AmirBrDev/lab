import matplotlib.pyplot as plt

DISTANCE_FROM_LIGATION_POINT = 10
DELIMITER = "\t"
SCORES_START_FIELD = 2
SEQ_FIELD = 1

negatives_scores = [0.0 for i in range(DISTANCE_FROM_LIGATION_POINT)]
positives_scores = [0.0 for i in range(DISTANCE_FROM_LIGATION_POINT)]

negative_read_count = 0

# calculate sums for each negative position
with open("/tmp/log_negative", "rb") as fl:

	for line in fl:
		sline = line.strip().split(DELIMITER)

		if len(sline[SEQ_FIELD]) < DISTANCE_FROM_LIGATION_POINT:
			continue

		negative_read_count += 1

		# Go over the last positions in the sequence
		for i in range(-1, - DISTANCE_FROM_LIGATION_POINT - 1, -1):
			negatives_scores[DISTANCE_FROM_LIGATION_POINT + i] += float(sline[i])


negatives_avg = [negatives_scores[i] / negative_read_count for i in range(DISTANCE_FROM_LIGATION_POINT)]

positive_read_count = 0

# calculate sums for each positive position
with open("/tmp/log_positive", "rb") as fl:

	for line in fl:
		sline = line.strip().split(DELIMITER)

		if len(sline[SEQ_FIELD]) < DISTANCE_FROM_LIGATION_POINT:
			continue

		positive_read_count += 1

		# Go over the last positions in the sequence
		for i in range(DISTANCE_FROM_LIGATION_POINT):
			positives_scores[i] += float(sline[i + SCORES_START_FIELD])

positives_avg = [positives_scores[i] / positive_read_count for i in range(DISTANCE_FROM_LIGATION_POINT)]

# Generate average base score plot
points = [[- DISTANCE_FROM_LIGATION_POINT + i, negatives_avg[i]] for i in range (DISTANCE_FROM_LIGATION_POINT)] + \
		 [[i + 1, positives_avg[i]] for i in range (DISTANCE_FROM_LIGATION_POINT)]

print points
plt.rcParams.update({"font.size": 18})
plt.xlim(-DISTANCE_FROM_LIGATION_POINT - 1, DISTANCE_FROM_LIGATION_POINT + 1)
plt.ylim(0, 1)
plt.hold(True)

for pt in points:
	plt.plot([pt[0], pt[0]], [0,pt[1]], "b")

plt.savefig("/tmp/random_ends_scores.png")