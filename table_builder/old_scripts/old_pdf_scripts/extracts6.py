

def read_to(fl, name):
    line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')
    while (not line.startswith(name)):

        line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')


fl = open("/Users/amirbar/Desktop/s6.txt", "rb")

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
    if i % 2 == 0:
        s5.append(lst)
        lst = []

read_to(fl, "Orientationb")

i = 0
line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')

while (not line.startswith("EOF")):

    if len(line) == 0:
        line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')
        continue

    s5[i].append(line)
    line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')
    i += 1

start_index = len(s5)

read_to(fl, "start_end")

i = 0
lst = []
line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')

while (not line.startswith("EOF")):

    if len(line) == 0:
        line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')
        continue

    lst.append(line)

    line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')


    i += 1
    if i % 2 == 0:
        s5.append(lst)
        lst = []

read_to(fl, "direction")

i = start_index
line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')

while (not line.startswith("EOF")):

    if len(line) == 0:
        line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')
        continue

    s5[i].append(line)
    line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')
    i += 1

start_index = len(s5)

read_to(fl, "start")
lst = []
line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')

while (not line.startswith("EOF")):

    if len(line) == 0:
        line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')
        continue

    s5.append([line])

    line = fl.readline().replace('\xe2\x80\xa9', '').replace('\n', '')


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

read_to(fl, "direction")

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


fl = open("/Users/amirbar/Desktop/s6_parsed.txt", "wb")

for lst in s5:
    fl.write("%s\n" % lst)

fl.close()