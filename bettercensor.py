#! /usr/bin/python

#Takes a string as input and replaces all incidents of "bad words" with "good words", as defined by a dictionary built-in to this program

def bettercensor(text):
    counterparts = {'ass': 'butt', 'shit': 'shoot', 'motherfucking': 'monkey fighting', 'bitch': 'witch'}
    sentence = text.split()
    newsentlist = []
    for term in range(len(sentence)):
        if sentence[term] in counterparts:
            newsentlist.append(counterparts[sentence[term]])
        else:
            newsentlist.append(sentence[term])
    newsent = ''
    for term in range(len(sentence)):
        if term == len(sentence) - 1:
            newsent+= newsentlist[term]
        else:
            newsent+= newsentlist[term] + ' '
    return newsent


bettercensor()
