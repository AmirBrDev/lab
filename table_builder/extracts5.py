

def read_to(fl, name):
    line = fl.readline().replace('\xe2\x80\xa9\n', '')
    while (not line.startswith(name)):

        line = fl.readline().replace('\xe2\x80\xa9\n', '')


fl = open("/Users/amirbar/Desktop/s5.txt", "rb")

line  = fl.readline().replace('\xe2\x80\xa9\n', '')

s5 = []
i = 0
lst = []

while (not line.startswith("EOF")):

    if len(line) == 0:
        line = fl.readline().replace('\xe2\x80\xa9\n', '').replace('\n', '')
        continue

    lst.append(line)

    line = fl.readline().replace('\xe2\x80\xa9\n', '')


    i += 1
    if i % 3 == 0:
        s5.append(lst)
        lst = []


read_to(fl, "sRNA")

i = 0
line = fl.readline().replace('\xe2\x80\xa9\n', '').replace('\n', '')

while (not line.startswith("EOF")):

    if len(line) == 0:
        line = fl.readline().replace('\xe2\x80\xa9\n', '').replace('\n', '')
        continue

    s5[i].append(line)
    line = fl.readline().replace('\xe2\x80\xa9\n', '').replace('\n', '')
    i += 1

read_to(fl, "MEV")

i = 0
line = fl.readline().replace('\xe2\x80\xa9\n', '').replace('\n', '')

while (not line.startswith("EOF")):

    if len(line) == 0:
        line = fl.readline().replace('\xe2\x80\xa9\n', '').replace('\n', '')
        continue

    s5[i].append(line)
    line = fl.readline().replace('\xe2\x80\xa9\n', '').replace('\n', '')
    i += 1

read_to(fl, "PRMc")

i = 0
line = fl.readline().replace('\xe2\x80\xa9\n', '').replace('\n', '')

while (not line.startswith("EOF")):

    if len(line) == 0:
        line = fl.readline().replace('\xe2\x80\xa9\n', '').replace('\n', '')
        continue

    s5[i].append(line)
    line = fl.readline().replace('\xe2\x80\xa9\n', '').replace('\n', '')
    i += 1

start_index = len(s5)

read_to(fl, "start")

line = fl.readline().replace('\xe2\x80\xa9\n', '').replace('\n', '')

while (not line.startswith("EOF")):

    if len(line) == 0:
        line = fl.readline().replace('\xe2\x80\xa9\n', '').replace('\n', '')
        continue

    s5.append([line])
    line = fl.readline().replace('\xe2\x80\xa9\n', '').replace('\n', '')

read_to(fl, "end")

line = fl.readline().replace('\xe2\x80\xa9\n', '').replace('\n', '')

i = start_index
while (not line.startswith("EOF")):

    if len(line) == 0:
        line = fl.readline().replace('\xe2\x80\xa9\n', '').replace('\n', '')
        continue

    s5[i].append(line)
    line = fl.readline().replace('\xe2\x80\xa9\n', '').replace('\n', '')
    i += 1

read_to(fl, "length")

line = fl.readline().replace('\xe2\x80\xa9\n', '').replace('\n', '')

i = start_index
while (not line.startswith("EOF")):

    if len(line) == 0:
        line = fl.readline().replace('\xe2\x80\xa9\n', '').replace('\n', '')
        continue

    s5[i].append(line)
    line = fl.readline().replace('\xe2\x80\xa9\n', '').replace('\n', '')
    i += 1

read_to(fl, "sRNA")

line = fl.readline().replace('\xe2\x80\xa9\n', '').replace('\n', '')

i = start_index
while (not line.startswith("EOF")):

    if len(line) == 0:
        line = fl.readline().replace('\xe2\x80\xa9\n', '').replace('\n', '')
        continue

    s5[i].append(line)
    line = fl.readline().replace('\xe2\x80\xa9\n', '').replace('\n', '')
    i += 1

read_to(fl, "MEV")

line = fl.readline().replace('\xe2\x80\xa9\n', '').replace('\n', '')

i = start_index
while (not line.startswith("EOF")):

    if len(line) == 0:
        line = fl.readline().replace('\xe2\x80\xa9\n', '').replace('\n', '')
        continue

    s5[i].append(line)
    line = fl.readline().replace('\xe2\x80\xa9\n', '').replace('\n', '')
    i += 1


read_to(fl, "PRM")

line = fl.readline().replace('\xe2\x80\xa9\n', '').replace('\n', '')

i = start_index
while (not line.startswith("EOF")):

    if len(line) == 0:
        line = fl.readline().replace('\xe2\x80\xa9\n', '').replace('\n', '')
        continue

    s5[i].append(line)
    line = fl.readline().replace('\xe2\x80\xa9\n', '').replace('\n', '')
    i += 1

fl.close()
print "finished"

for lst in s5:
    print lst


fl = open("/Users/amirbar/Desktop/s5_parsed.txt", "wb")

for lst in s5:
    fl.write("%s\n" % lst)

fl.close()