from Table import Table
from random import randint
import numpy as np

def random_names_for_training(pool, training_fraction):
	result = []

	tmp = [name for name in pool]

	count  = int(len(pool) * training_fraction)

	print count
	while len(result) < count:
		#print len(result), len(tmp)
		name = tmp[randint(0, len(tmp) - 1)]
		result.append(name)
		tmp.remove(name)

	return result

def retrieve_data(data_pool, subset):
	vectors = []

	names_to_rows = {row["name"]: row for row in data_pool}

	for name in subset:

		if name in names_to_rows.keys():

			row = names_to_rows[name]

			present_count = 0.
			second_precentage = 0.
			il_percentage = row["signif_chimeras_of_iron_limitation_cl.as_rna2_percentage"]

			if il_percentage != "-":
				present_count += 1
				second_precentage += float(il_percentage)

			log_percentage = row["signif_chimeras_of_log_phase_cl.as_rna2_percentage"]

			if log_percentage != "-":
				present_count += 1
				second_precentage += float(log_percentage)

			stat_percentage = row["signif_chimeras_of_stationary_cl.as_rna2_percentage"]

			if stat_percentage != "-":
				present_count += 1
				second_precentage += float(stat_percentage)

			vector = np.array([float(row["total_interactions"]),
							   float(row["tb_srna_targets"]),
							   float(row["mrna_5utr_targets"]),
							   float(row["max_poly_u_length"]),
							   second_precentage / present_count])

			vector = np.array([float(vector[i]) for i in range(vector.shape[0])])

			vectors.append(vector)
		else:
			print name

	return vectors

def load_file(file_path):
	data_pool = Table(file_path)	

	names = [row["name"] for row in data_pool]
	other_fields = [row for row in data_pool]

	return names, retrieve_data(data_pool, names), other_fields


def retrieve_training_set():

	training_pool = Table("training_pool_3_targets.table")

	tmp = []
	for name in [row["name"] for row in training_pool]:
		if name in tmp:
			print name
		else:
			tmp.append(name)

	positive_training_pool = []
	negative_training_pool = []

	training_fraction = 0.7

	print "<---", "cpxR" in [row["name"] for row in training_pool]

	for row in training_pool:
		if int(row["label"]) < 0:
			negative_training_pool.append(row["name"])
		else:
			positive_training_pool.append(row["name"])

	data_pool = Table("table_s6_range_0.csv")

	positive_training_names = \
		random_names_for_training(positive_training_pool, 
								  training_fraction)

	negative_training_names = \
		random_names_for_training(negative_training_pool, 
								  training_fraction)

	positive_validation_names = \
		[name for name in positive_training_pool if name not in positive_training_names]
	negative_validation_names = \
		[name for name in negative_training_pool if name not in negative_training_names]

	training_set = retrieve_data(data_pool, positive_training_names)
	training_set.extend(retrieve_data(data_pool, negative_training_names))
	training_labels = [1] * len(positive_training_names)
	training_labels.extend([0] * len(negative_training_names))

	validation_set = retrieve_data(data_pool, positive_validation_names)
	validation_set.extend(retrieve_data(data_pool, negative_validation_names))
	validation_labels = [1] * len(positive_validation_names)
	validation_labels.extend([0] * len(negative_validation_names))

	training_set = np.array(training_set)
	training_labels = np.array(training_labels)
	validation_set = np.array(validation_set)
	validation_labels = np.array(validation_labels)

	training_names = positive_training_names + negative_training_names
	validation_names = positive_validation_names + negative_validation_names

	print len(training_set), len(training_labels)

	return training_set, training_labels, validation_set, validation_labels, training_names, validation_names
	
if __name__ == "__main__":

	print retrieve_training_set()

