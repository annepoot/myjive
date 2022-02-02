def parse(fname):
    global i
    global sp

    with open(fname, 'r') as input:
        filestr = input.read().replace('\n', '').replace('\t', '')

    data = {}

    i = 0
    sp = filestr.split(';')

    while i < len(sp):
        line = sp[i].replace(' ','')
        if '{' in line:
            key = line.split('={')[0]
            newline = '={'.join(line.split('={')[1:])
            data[key] = readLevel(newline)
        elif '=' in line:
            [key, value] = line.split('=')
            subdata[key] = value
        elif line != '':
            raise RuntimeError ('Unable to parse: %s' % line)

        i = i + 1

    return data


def readLevel(line):
    global i
    global sp
    subdata = {}

    while True:
        if '{' in line: 
            key = line.split('={')[0]
            newline = '={'.join(line.split('={')[1:])
            subdata[key] = readLevel(newline)
        elif '}' in line:
            return subdata
        elif '=' in line:
            [key, value] = line.split('=')
            subdata[key] = value
        elif line != '':
            raise RuntimeError ('Unable to parse: %s' % line)

        i = i + 1

        if i == len(sp):
            raise RuntimeError('EOF reached while parsing an input block. Did you forget to close a bracket?')

        line = sp[i].replace(' ','')

