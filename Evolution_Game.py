#! /usr/bin/python
import os
import random
def Evolution_Game():
    X = [] #"Females" of species 1
    Y = [] #"Males" of species 1

    W = []
    Z = []

    X_0 = []
    Y_0 = []

    W_0 = []
    Z_0 = []

    base = 17.
    sigma = 4.
    X_gap = 4. #How far below average an X-type must be to be unacceptable
    Y_gap = 2. #How far below average a Y-type must be to be unacceptable

    W_gap = 5.
    Z_gap = 7.


    for x in range(0, 100):
        X_0.append(random.gauss(base,sigma))
        Y_0.append(random.gauss(base,sigma))
        W_0.append(random.gauss(base,sigma))
        Z_0.append(random.gauss(base,sigma))

    X.append(X_0)
    Y.append(Y_0)
    W.append(W_0)
    Z.append(Z_0)


    x_min = base - X_gap
    y_min = base - Y_gap

    X_total = [sum(X[0])]
    Y_total = [sum(Y[0])]

    x_averages = [base]
    y_averages = [base]

    x_mins = [x_min]
    y_mins = [y_min]

    pop_x = [len(X[0])]
    pop_y = [len(Y[0])]

    print "This is the initial X population %f" % pop_x[0]

    #X_1 = []
    #Y_1 = []

    for i in range(1,101):
        tmp_x = []
        tmp_y = []
        
        if pop_x[i-1] < pop_y[i-1]:
            size = pop_x[i-1]
        else:
            size = pop_y[i-1]

        for j in xrange(size):  
            if (Y[i-1][j] > y_mins[i-1]) & (X[i-1][j] > x_mins[i-1]):
                offcount = random.randint(4,8) #This determines the reproductive rate of the population
                for z in xrange(offcount):
                    bio = random.randint(1,2)
                    if bio == 1:
                        bar_x = X[i-1][j]*0.9 + Y[i-1][j]*0.1
                        new_x = random.gauss(bar_x, 0.1*bar_x) #Should be updated to have the actual sigma of the parents
                        tmp_x.append(new_x)
                    else:
                        bar_y = X[i-1][j]*0.1 + Y[i-1][j]*0.9
                        new_y = random.gauss(bar_y, 0.1*bar_y) #See above comment
                        tmp_y.append(new_y)
                        
        X.append(tmp_x)
        Y.append(tmp_y)
        
        pop_x.append(len(tmp_x))
        pop_y.append(len(tmp_y))
                        
        #print "This is the X population %f" % pop_x[i]
        #print "This is the Y population %f" % pop_y[i]

        X_total.append(sum(tmp_x))
        Y_total.append(sum(tmp_y))

        #print "This is tmp_x ", tmp_x
        #print "This is tmp_y ", tmp_y

        if (pop_x[i] > 0) & (pop_y[i] > 0):
            x_averages.append(X_total[i]/float(pop_x[i]))
            y_averages.append(Y_total[i]/float(pop_y[i]))
            x_mins.append(x_averages[i] - X_gap)
            y_mins.append(y_averages[i] - Y_gap)
        else:
            x_averages.append(0)
            y_averages.append(0)
            x_mins.append(x_averages[i-1] - X_gap)
            y_mins.append(y_averages[i-1] - Y_gap)

        if (i % 10) == 0:
            print "This is generation %f" % i
            print "The current X population is %f" % pop_x[i]
            print "The current Y population is %f" % pop_y[i]
            print "The current average x is %f" % x_averages[i]
            print "The current average y is %f" % y_averages[i]



    #for j in xrange(len(X_0)):  
    #    if (Y_0[j] > y_min) & (X_0[j] > x_min):
    #        offcount = random.randint(1,4)
    #        for z in range(offcount):
    #            bio = random.randint(1,2)
    #            if bio == 1:
    #                new_x = X_0[j]*0.9 + Y_0[j]*0.1
    #                X_1.append(new_x)
    #            else:
    #                new_y = X_0[j]*0.1 + Y_0[j]*0.9
    #                Y_1.append(new_y)

    #X_total = sum(X_1)
    #Y_total = sum(Y_1)
    #X1_mean = X_total/len(X_1)
    #Y1_mean = Y_total/len(Y_1)

    #print "The new X-mean is %f" % X1_mean
    #print "The new Y-mean is %f" % Y1_mean
    
Evolution_Game()
