def parse(fname):
    global i
    global sp

    with open(fname, 'r') as input:
        filestr = input.read().replace('\n', '').replace('\t', '')

    data = {}

    i = 0
    sp = filestr.split(';')

    while i < len(sp):
        line = sp[i]
        i = i + 1
        if '{' in line.replace(' ', ''):
            [key, value] = line.replace(' ', '').split('={')
            data[key] = readLevel()
        elif 'include' in line:
            data[line] = None

    return data


def readLevel():
    global i
    global sp
    subdata = {}

    if 'include' in sp[i - 1]:
        subdata[sp[i - 1].split('{')[1]] = None
    elif len(sp[i - 1].split('=')) == 2:
        [key, value] = sp[i - 1].replace(' ', '').split('={')
        subdata[key] = value
    else:
        [key, value] = sp[i - 1].replace(' ', '').split('{')[1].split('=')
        subdata[key.replace(' ', '')] = value

    while True:
        line = sp[i]
        i = i + 1
        if '{' in line.replace(' ', ''):
            [key, value] = line.replace(' ', '').split('={')
            subdata[key] = readLevel()
        elif '}' in line:
            return subdata
        elif 'include' in line:
            subdata[line] = None
        else:
            [key, value] = line.replace(' ', '').split('=')
            subdata[key.replace(' ', '')] = value

