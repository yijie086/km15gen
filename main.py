from km15gen import genOneEvent
import argparse
import numpy as np
import time
import os
import subprocess

M  = 0.938272081
x1 = 1/2/M/8.604
x2 = 1/(5-M**2)
x3 = (10.604/8.604-1)/M*10.604* (1-np.cos(np.radians(35)))
x4 = (1-(4-M**2)/2/10.604/M)/(1+(4-M**2)/2/10.604**2/(1-np.cos(np.radians(35))))

y1 = 1
y2 = 1.456
y3 = 2.510
y4 = 4.326
y5 = 7.671

c0 = y2/2/M/8.604
d0 = 1/(1+(4-M*M)/y2)
c1 = np.sqrt(y2*y3)/2/M/8.604
d1 = 1/(1+(4-M*M)/np.sqrt(y2*y3))
c2 = y3/2/M/8.604
d2 =  1/(1+(4-M*M)/y3)
c3 = np.sqrt(y3*y4)/2/M/8.604
d3 = 1/(1+(4-M*M)/np.sqrt(y3*y4))
c4 = y4/2/M/8.604
d4 = 1/(1+(4-M*M)/y4)
c5 = np.sqrt(y4*y5)/2/M/8.604
d5 = 1/(1+(4-M*M)/np.sqrt(y4*y5))
c6 = y5/2/M/8.604
d6 = 1/(1+(4-M*M)/y5)

newxBbins = [x1, c0, c1, c2, c3, c4, d2]
newQ2bins = [y1, y2, np.sqrt(y2*y3), y3, np.sqrt(y3*y4), y4, np.sqrt(y4*y5)]
newxBbins2 = [x1, c0, c1, c2, c3, c4,                             c5, d2, d4]
newQ2bins2 = [y1, 1.2, y2, np.sqrt(y2*y3), y3, np.sqrt(y3*y4)    , y4, np.sqrt(y4*y5), 7]
newtbins = [0.11, 0.15, 0.25, 0.4, 0.6, 0.8, 1.0, 1.25, 1.5, 1.79]

def main(args):

  Ed    = args.Ed
  if args.bin:
    bin_scheme = np.loadtxt("/work/clas12/sangbaek/km15gen/bin_scheme.csv")
    xBmin, xBmax, Q2min, Q2max, tmin, tmax = bin_scheme[args.bin - 1]
  else:
    xBmin = args.xBmin
    xBmax = args.xBmax
    Q2min = args.Q2min
    Q2max = args.Q2max
    tmin  = args.tmin
    tmax  = args.tmax
  ymin  = args.ymin
  ymax  = args.ymax
  w2min  = args.w2min
  rad = int(args.radgen)
  trig = args.trig
  filename = args.fname

  now = time.time()
  with open("{}.dat".format(filename), "w") as file_out:
    file_out.write("")

  if args.model == 'km15':
    num  = 0
    while num <trig:
      result   = genOneEvent(xBmin, xBmax, Q2min, Q2max, tmin, tmax, ymin, ymax, w2min, 0,  rad = rad, Ed = Ed, filename = filename)
      num      = num + result

    later = time.time()
    print("The time spent in generating events: {:.3f} s".format(later-now))
    now = time.time()

    # remove the last line break
    file = open("{}.dat".format(filename), "r+")
    lines = file.readlines()
    lines[-1] = lines[-1][:-1]
    file = open("{}.dat".format(filename), "w")
    file.writelines(lines)

  elif args.model == 'bh':
    dvcsgen_commands = ["/work/clas12/sangbaek/dvcsgen/dvcsgen", "--docker", "--trig", "{}".format(trig), "--beam", "{:.3f}".format(Ed),
      "--x", "{:.3f}".format(xBmin), "{:.3f}".format(xBmax),
      "--q2", "{:.3f}".format(Q2min), "{:.3f}".format(Q2max),
      "--t", "{:.3f}".format(tmin), "{:.3f}".format(tmax),
      "--gpd", "101", "--y", "{:.3f}".format(ymin), "{:.3f}".format(ymax), "--w", "{:.3f}".format(w2min),
      "--raster", "0.025", "--writef", "2", "--globalfit", "--ycol", "0.0005"]
    if rad:
      dvcsgen_commands.extend(["--radgen", "--vv2cut", "0.6", "--delta", "0.1", "--radstable"])
    dvcsgen_commands.extend(["--bh", "1"])
    subprocess.run(dvcsgen_commands)

  elif args.model == 'vgg':
    dvcsgen_commands = ["/work/clas12/sangbaek/dvcsgen/dvcsgen", "--docker", "--trig", "{}".format(trig), "--beam", "{:.3f}".format(Ed),
      "--x", "{:.3f}".format(xBmin), "{:.3f}".format(xBmax),
      "--q2", "{:.3f}".format(Q2min), "{:.3f}".format(Q2max),
      "--t", "{:.3f}".format(tmin), "{:.3f}".format(tmax),
      "--gpd", "101", "--y", "{:.3f}".format(ymin), "{:.3f}".format(ymax), "--w", "{:.3f}".format(w2min),
      "--raster", "0.025", "--writef", "2", "--globalfit", "--ycol", "0.0005"]
    if rad:
      dvcsgen_commands.extend(["--radgen", "--vv2cut", "0.6", "--delta", "0.1", "--radstable"])
    dvcsgen_commands.extend(["--bh", "3"])
    subprocess.run(dvcsgen_commands)

  elif args.model == 'pi0':

    aao_gen_commands = ['/work/clas12/sangbaek/aao_gen/gen_wrapper/src/aao_gen.py', '--generator_type', 'rad', 
    '--input_filename_rad', './aao_rad_input.inp', '--input_filename_norad', './aao_norad_input.inp',
    '--flag_ehel', '1', '--ebeam', '{:.3f}'.format(Ed), '--q2min', '{:.3f}'.format(Q2min), '--q2max', '{:.3f}'.format(Q2max),
    '--epmin', '0.1', '--epmax', '10.604', '--fmcall', '1.0', '--boso', '1',
    '--seed', '0', '--trig', '{}'.format(trig), '--epirea', '1', '--physics_model_rad', '5',
    '--int_region', '.20 .12 .20 .20', '--npart_rad', '4', '--sigr_max_mult', '0.0',
    '--sigr_max', '0.005', '--model_5_min_W', '3.5721', '--rad_emin', '0.005', '--err_max', '0.2',
    '--target_len', '5', '--target_rad', '0.43', '--cord_x', '0.0', '--cord_y', '0.0', '--cord_z', '-3',
    '--physics_model_norad', '5', '--npart_norad', '3', '--input_exe_path', '/work/clas12/sangbaek/aao_gen/gen_wrapper/src/aao_input_file_maker.py',
    '--precision', '5', '--maxloops', '10', '--generator_exe_path', '/work/clas12/sangbaek/aao_gen/aao_rad/build/aao_rad',
    '--xBmin', '{:.3f}'.format(xBmin), '--xBmax', '{:.3f}'.format(xBmax), '--w2min', '{:.3f}'.format(w2min),
    '--w2max', '50.0', '--tmin', '{:.3f}'.format(tmin), '--tmax', '{:.3f}'.format(tmax), '--filter_infile',
    './aao_rad.lund', '--filter_outfile', './aao_gen.dat', '--filter_exe_path', '/work/clas12/sangbaek/aao_gen/gen_wrapper/src/lund_filter.py',
    '--outdir', '.']
    if rad:
      pass
    else:
      print("non-radiative generator for pi0 is not needed. Exit.")
      exit()
    subprocess.run(aao_gen_commands)

  # if args.find_max:
  #   xs_array = []
  #   while len(xs_array) < args.trig:
  #     xs = bhdvcs(xBmin, xBmax, Q2min, Q2max, tmin, tmax, ymin, ymax, xs_only = 1)
  #     if xs:
  #       xs_array.append(xs)

  #   print(np.max(xs_array))
if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument("-Ed", "--Ed", type = float, default = 10.604)
  parser.add_argument("-trig", "--trig", type = int, default = 1)
  parser.add_argument("-fname", "--fname", type = str, default = "km15gen")
  parser.add_argument("-bin", "--bin", type = int, default = 0)
  parser.add_argument("-model", "--model", type = str, default = 'km15')
  parser.add_argument("-xBmin", "--xBmin", type = float, default = 0.05)
  parser.add_argument("-xBmax", "--xBmax", type = float, default = 0.75)
  parser.add_argument("-Q2min", "--Q2min", type = float, default = 0.9)
  parser.add_argument("-Q2max", "--Q2max", type = float, default = 11)
  parser.add_argument("-tmin", "--tmin", type = float, default = 0.085)
  parser.add_argument("-tmax", "--tmax", type = float, default = 1.79)
  parser.add_argument("-ymin", "--ymin", type = float, default = 0.19)
  parser.add_argument("-ymax", "--ymax", type = float, default = 0.85)
  parser.add_argument("-w2min", "--w2min", type = float, default = 3.61)
  parser.add_argument("-radgen", "--radgen", action = 'store_true')
  args = parser.parse_args()

  main(args)
