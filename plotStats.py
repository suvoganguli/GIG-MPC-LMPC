import matplotlib.pyplot as plt
import numpy as np

# delim="\t"
def getColumns(inFile, delim=" ", header=True):
    """
    Get columns of data from inFile. The order of the rows is respected

    :param inFile: column file separated by delim
    :param header: if True the first line will be considered a header line
    :returns: a tuple of 2 dicts (cols, indexToName). cols dict has keys that
    are headings in the inFile, and values are a list of all the entries in that
    column. indexToName dict maps column index to names that are used as keys in
    the cols dict. The names are the same as the headings used in inFile. If
    header is False, then column indices (starting from 0) are used for the
    heading names (i.e. the keys in the cols dict)
    """
    cols = {}
    indexToName = {}
    for lineNum, line in enumerate(inFile):
        if lineNum == 0:
            headings = line.split(delim)
            i = 0
            for heading in headings:
                heading = heading.strip()
                if header:
                    cols[heading] = []
                    indexToName[i] = heading
                else:
                    # in this case the heading is actually just a cell
                    cols[i] = [heading]
                    indexToName[i] = i
                i += 1
        else:
            cells = line.split(delim)
            i = 0
            for cell in cells:
                cell = cell.strip()
                cols[indexToName[i]] += [cell]
                i += 1

    return cols, indexToName

filenames = ['logFile_nstates4_N6.txt','logFile_nstates4_N8.txt',
             'logFile_nstates4_N10.txt', 'logFile_nstates4_N12.txt']

n = len(filenames)
cputime_vec = np.zeros(n)
N_vec = [6,8,10,12]


plt.figure(1)
for k in range(n):
    f = file(filenames[k], 'r')
    cols, indexToName = getColumns(f, header=False)
    f.close()

    iter = np.array(cols[0]).astype(np.int)
    cputime = np.array(cols[1]).astype(np.float)

    cputime_vec[k] = np.mean(cputime)

    plt.plot(iter,cputime)
    plt.grid('True')
    plt.xlabel('Iterations')
    plt.ylabel('CPU time [sec]')
    plt.title('nstates = 4')
    # plt.gca().set_adjustable('box')
    # ax = plt.gca()
    # ax.set_ylim(0,5))
plt.legend(['N = 6','N = 8','N = 10','N = 12'])

plt.figure(2)
plt.plot(N_vec,cputime_vec, marker = 'x')
plt.xlabel('N')
plt.ylabel('Average CPU time [sec]')
plt.grid('True')

plt.show()


