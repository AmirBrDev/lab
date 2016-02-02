def sigmoid(X):
	return 1./ (1. + numpy.exp(- X))


def cost(theta, X, y):

	p_1 = sigmoid(numpy.dot(X, theta)) # predicted probability of label 1

	log_l = []

	for i, p_i in enumerate(p_1):

		if y[i] == 1:
			log_l.append(numpy.log(p_i))

		if y[i] == 0:
			log_l.append(numpy.log(1-p_i))

	log_l = numpy.array(log_l)

	return - log_l.mean()

def grad(theta, X, y):
	p_1 = sigmoid(numpy.dot(X, theta))
	error = p_1 - y # difference between label and prediction
	grad = numpy.dot(error, X_1) / y.size # gradient vector

	return grad

def predict(theta, X):
	p_1 = sigmoid(numpy.dot(X, theta))
	return p_1 > 0.5

def odds(theta, X):
	return sigmoid(numpy.dot(X, theta))

import numpy
import scipy.optimize as opt


def test():
	# prefix an extra column of ones to the feature matrix (for intercept term)
	#theta = 0.1* numpy.random.randn(3)
	theta = numpy.array([1., 1., 1.])

	training_set = [([1.,1.], 0.),
					([10.,20.], 0.),
					([49., 49.], 0.),
					([49., 50.], 0.),
					([50., 49.], 0.),
					([50., 51.], 1.),
					([51., 50.], 1.),
					([51.,51.], 1.),
					([100.,100.], 1.),
					([0., 100.], 0.),
					([100., 0.], 0.),
					([51., 100.], 1.),
					([100., 51.], 1.),
					([50., 60.], 1.)]

	X = numpy.array([x for x, y in training_set])
	y = numpy.array([y for x, y in training_set])

	X_norm = (X - X.mean(axis=0)) / X.max(axis=0)

	print X_norm

	X_1 = numpy.append( numpy.ones((X_norm.shape[0], 1)), X_norm, axis=1)

	theta_1 = opt.fmin_bfgs(cost, theta, fprime=grad, args=(X_1, y))

	test_set = [[50.5, 50.5]]
	test_set = (numpy.array(test_set) - X.mean(axis=0)) / X.max(axis=0)
	print test_set

	test_set_1 = numpy.append( numpy.ones((test_set.shape[0], 1)), test_set, axis=1)

	for i in test_set_1:
		print i, predict(numpy.array(theta_1), i), odds(numpy.array(theta_1), i)

from training_loader import retrieve_training_set, load_file

print "*" * 100
print """Running logistic regression"""
print "*" * 100

X, y, val_x, val_y, training_names, validation_names= retrieve_training_set()

X_norm = (X - X.mean(axis=0)) / X.max(axis=0)

X_1 = numpy.append( numpy.ones((X_norm.shape[0], 1)), X_norm, axis=1)

theta = numpy.array([1.] * (X.shape[1] + 1))
theta_1 = opt.fmin_bfgs(cost, theta, fprime=grad, args=(X_1, y))

val_x_norm = (val_x - X.mean(axis=0)) / X.max(axis=0)
val_x_1 = numpy.append( numpy.ones((val_x_norm.shape[0], 1)), val_x_norm, axis=1)

print "#" * 80
print "Running validation"

for i, name in zip(val_x_1, validation_names):
	print name, predict(numpy.array(theta_1), i), odds(numpy.array(theta_1), i)

print "theta:", theta_1
print "Validation Finished"

print "#" * 80
print "Analyze data:"
names, vectors, other_fields = load_file("./table_s6_range_0_above3_targets.csv")

vectors_norm = (vectors - X.mean(axis=0)) / X.max(axis=0)
vectors_padded = numpy.append( numpy.ones((vectors_norm.shape[0], 1)), vectors_norm, axis=1)

with open("results.txt", "wb") as fl:

	for name, vector, other in zip(names, vectors_padded, other_fields):
		fl.write("%s\t%s\t%s\t%s\t%s\n" % \
			(name, 
			 predict(numpy.array(theta_1), vector), 
			 odds(numpy.array(theta_1), vector),
			 other["total_interactions"],
			 other["mrna_5utr_targets"]))

	fl.close()

print "*" * 100
