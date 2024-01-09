import sys
import numpy as np
import matplotlib.pyplot as plt


def read_data (input_file):
	data = np.genfromtxt(input_file)
	return data[:,0], data[:,1]

def plot_data(x100, y100, x200, y200, x400, y400):
    plt.grid()
    plt.plot(x100, y100, label="100x100 (h=0.5)", c="red", linewidth=1.0)
    plt.plot(x200, y200, label="200x200 (h=0.25)", c="green", linewidth=1.0)
    plt.plot(x400, y400, label="400x400 (h=0.125)", c="blue", linewidth=1.0)
    plt.ylabel("Ce", fontsize=15)
    plt.xlabel("time (s)", fontsize=15)
    plt.title("Mean extravascular concentration overtime", fontsize=14)
    plt.legend(loc=0,fontsize=10)
    #plt.show()
    #plt.savefig("extra_vascular_concentration.pdf")
    plt.savefig("extra_vascular_concentration.png", dpi=300)

def main():
    if len(sys.argv) != 4:
        print("-----------------------------------------------------------------------------------------------------------")
        print("Usage:> python %s <input_Ce_vent100x100> <input_Ce_vent50x50> <input_Ce_vent25x25>" % sys.argv[0])
        print("-----------------------------------------------------------------------------------------------------------")
        print("<input_Ce_vent100x100> = Input file with the extravascular concentration for the 100x100 grid")
        print("<input_Ce_vent200x200> = Input file with the extravascular concentration for the 200x200 grid")
        print("<input_Ce_vent400x400> = Input file with the extravascular concentration for the 400x400 grid")
        print("-----------------------------------------------------------------------------------------------------------")
        return 1

    input_Ce_vent100x100_file = sys.argv[1]
    input_Ce_vent200x200_file = sys.argv[2]
    input_Ce_vent400x400_file = sys.argv[3]

    x100, y100 = read_data(input_Ce_vent100x100_file)
    x200, y200 = read_data(input_Ce_vent200x200_file)
    x400, y400 = read_data(input_Ce_vent400x400_file)
    plot_data(x100, y100, x200, y200, x400, y400)

if __name__ == "__main__":
	main()
