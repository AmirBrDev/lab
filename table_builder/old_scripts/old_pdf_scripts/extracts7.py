__author__ = 'amirbar'



def read_to(fl, name):
    line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')
    while (not line.startswith(name)):

        line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')


fl = open("/Users/amirbar/Desktop/s7.txt", "rb")

line  = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')

s5 = []
i = 0
lst = []

while (not line.startswith("EOF")):

    if len(line) == 0:
        line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')
        continue

    lst.append(line)

    line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')


    i += 1
    if i % 8 == 0:
        s5.append(lst)
        lst = []

start_index = len(s5)

read_to(fl, "name")
lst = []
line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')

while (not line.startswith("EOF")):

    if len(line) == 0:
        line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')
        continue

    s5.append([line])

    line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')


read_to(fl, "start")

i = start_index
line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')

while (not line.startswith("EOF")):

    if len(line) == 0:
        line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')
        continue

    s5[i].append(line)
    line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')
    i += 1

read_to(fl, "end")

i = start_index
line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')

while (not line.startswith("EOF")):

    if len(line) == 0:
        line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')
        continue

    s5[i].append(line)
    line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')
    i += 1

read_to(fl, "length")

i = start_index
line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')

while (not line.startswith("EOF")):

    if len(line) == 0:
        line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')
        continue

    s5[i].append(line)
    line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')
    i += 1


read_to(fl, "strand")

i = start_index
line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')

while (not line.startswith("EOF")):

    if len(line) == 0:
        line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')
        continue

    s5[i].append(line)
    line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')
    i += 1

read_to(fl, "MEV")

i = start_index
line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')

while (not line.startswith("EOF")):

    if len(line) == 0:
        line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')
        continue

    s5[i].append(line)
    line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')
    i += 1

read_to(fl, "PRM")

i = start_index
line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')

while (not line.startswith("EOF")):

    if len(line) == 0:
        line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')
        continue

    s5[i].append(line)
    line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')
    i += 1

read_to(fl, "source")

i = start_index
line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')

while (not line.startswith("EOF")):

    if len(line) == 0:
        line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')
        continue

    s5[i].append(line)
    line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')
    i += 1


fl.close()
print "finished"

for lst in s5:
    print lst


fl = open("/Users/amirbar/Desktop/s7_parsed.txt", "wb")

for lst in s5:
    fl.write("%s\n" % lst)

fl.close()