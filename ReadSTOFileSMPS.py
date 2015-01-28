import re

# Read the .sto file in SMPS format
def read_sto_file_smps(file_name,nScen,nClient,nServer):
    file = open(file_name, 'r')
    s = 1
    h = [[0 for col in range(nClient)] for row in range(nScen)]
    line = file.next()
    while(not(line.startswith('ENDATA'))):
        if line.strip().upper().startswith('SC SCEN'+str(s)):
            line = file.next()
            while(line.strip().startswith('RHS')):
                i = line.split()[1]
                i = int(i[1:]) - 1 - nServer - 1
                val = int(line.split()[2])
                #print s, i
                h[s-1][i] = val
                #print val
                #print line
                line = file.next()
            s += 1
        else:
            line = file.next()

    file.close()
    return h





