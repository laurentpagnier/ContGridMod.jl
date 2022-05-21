using HDF5

data = h5read("pantagruel.h5","/")
data["bus_coord"] = data["bus_coord"][:,2:-1:1]

h5write("panta.h5","gen_prim_ctrl", data["gen_prim_ctrl"])
h5write("panta.h5","load_freq_coef", data["load_freq_coef"])
h5write("panta.h5","branch", data["branch"])
h5write("panta.h5","bus_coord", data["bus_coord"])
h5write("panta.h5","gen", data["gen"])
h5write("panta.h5","gen_inertia", data["gen_inertia"])
h5write("panta.h5","bus", data["bus"])
