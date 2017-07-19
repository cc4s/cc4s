import numpy
import os

folder = os.path.dirname(__file__) or "."
print(folder)

file_name = folder+"/"+"SimlarityTransformedHamiltonianSD.dat"
data = numpy.loadtxt(file_name, skiprows=2)

plot = False
if plot:
    import matplotlib.pyplot as plt
    plt.matshow(data)
    fig = plt.gcf()
    plt.colorbar()
    plt.savefig(folder+"/matrix.pdf")
    plt.savefig("matrix.pdf")


# print("Determinant")
try:
    det = numpy.linalg.det(data)
except:
    det = 0

print("Determinant    %s" % det)

if det or True:
    eigva, eigve = numpy.linalg.eig(data)

    print("\n\nEigenvalues")
    for e in eigva:
        print(e)

    print("\n\nEigenvalues(sorted)")
    for j,e in enumerate(numpy.sort(eigva)):
        print("%s.  %s" % (j, e))

    # print("\n\nEigenvectors")
    # for ev in eigve:
        # print(ev)


#vim-run: python3 %
#vim-run: python % && mupdf matrix.pdf
