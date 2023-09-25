folder="../4x4x4x32/b2p44_new/"
measure="gf_afm/"
data=np.loadtxt(folder+measure+"Index_ov.txt")
ov_max, susy_max=np.loadtxt(folder+measure+"end_spectrum.txt")
plt.plot(data[1],data[0], label="$t=4.0$")
plt.scatter(susy_max[0],data[0,int(susy_max[1])], color="red", marker="v")

measure="gf_afm_2p0t/"
data=np.loadtxt(folder+measure+"Index_ov.txt")
ov_max, susy_max=np.loadtxt(folder+measure+"end_spectrum.txt")
plt.plot(data[1],data[0], label="$t=2.0$")
plt.scatter(susy_max[0],data[0,int(susy_max[1])], color="red", marker="v")

measure="gf_afm_1p125t/"
data=np.loadtxt(folder+measure+"Index_ov.txt")
ov_max, susy_max=np.loadtxt(folder+measure+"end_spectrum.txt")
plt.plot(data[1],data[0], label="$t=1.125$")
plt.scatter(susy_max[0],data[0,int(susy_max[1])], color="red", marker="v")

measure="gf_afm_0p5t/"
data=np.loadtxt(folder+measure+"Index_ov.txt")
ov_max, susy_max=np.loadtxt(folder+measure+"end_spectrum.txt")
plt.plot(data[1],data[0], label="$t=0.5$")
plt.scatter(susy_max[0],data[0,int(susy_max[1])], color="red", marker="v")

measure="gf_afm_0p25t/"
data=np.loadtxt(folder+measure+"Index_ov.txt")
ov_max, susy_max=np.loadtxt(folder+measure+"end_spectrum.txt")
plt.plot(data[1],data[0], label="$t=0.25$")
plt.scatter(susy_max[0],data[0,int(susy_max[1])], color="red", marker="v")

plt.legend(loc="upper right")
plt.xlabel(r'$\lambda_{max}$')
plt.title("Configurations with wrong index $\mathcal{O}$")
plt.savefig("./Index_ov.pdf")