#!/usr/bin/python

# numpy array capabilities are needed
import numpy as np

# os and glob for file operations
import os.path
import glob

# sys to read the arguments
import sys

# plotting library
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


##############################################################################
#
# Explain how this program works in general.
#
##############################################################################
general_man = 'This python script creates some plots for Naveeds\nmixing'\
+ ' layer problem. The syntax is\n\n  ./PlotResults <option>=<value> [DIR/FILE]'\
+ '\n\n'\
+ 'or\n\n  ./configure [-h, --help]\n\nto display this help.\n'\
+ '==================================================\n'\
+ 'Allowed <option>s and <value>s are:\n'

##############################################################################
#
# Explain how the directory or file is given.
#
##############################################################################
dirfile_man = '  [DIR/FILE]:\n    The last argument must be a directory or'\
+ ' a\n    single file. Relative paths are allowed, e.g.:\n\n'\
+ '      ./PlotResults --format=*_supg_*p2*.out .\n'

##############################################################################
#
# Explain how line style options work.
#
##############################################################################
styles_man = '  --styles:\n    To give the line styles for each file\n'\
+ '    use comma seperated values, as in matplotlib\n    plot() function,'\
+ ' e.g.:\n\n'\
+ '      --styles=g+,b-.,m--\n\n'
    
##############################################################################
#
# Explain how figures can be chosen.
#
##############################################################################
figures_man = '  --figures:\n    To specify which figures to plot,\n'\
+ '    use comma seperated values.\n    Default <value>s are 1,3,5,6.\n'\
+ '    Figure 1: y over u mean\n'\
+ '    Figure 2: x over v mean\n'\
+ '    Figure 3: y over u root mean squared\n'\
+ '    Figure 4: x over v root mean squared\n'\
+ '    Figure 5: kinetic energy over time\n'\
+ '    Figure 6: relative vorticity thickness over time\n'\
+ '    Example usage:\n\n'\
+ '      --figures=1,2,6\n'

##############################################################################
#
# Explain how file format works.
#
##############################################################################
format_man = '  --format:\n    To give a file format for the files\n'\
+ '    you want to plot.\n'\
+ '    This works as in the unix file system, e.g.:\n\n'\
+ '      --format=*_supg_*p2*.out\n'

##############################################################################
#
# THIS FUNCTION ALWAYS GOES FIRST IN THE MAIN PART.
#
# Check if the help is called or no arguments are given.
# This aborts everything else, just prints the help and writes no cache file.
#
##############################################################################
def check_for_help(opts):
  global general_man
  global dirfile_man
  global styles_man
  global figures_man
  global format_man
  
  begin_end_str = '==================================================\n'
  delim_str = '--------------------------------------------------\n'
  
  man_string = begin_end_str \
  + general_man \
  + dirfile_man + delim_str \
  + styles_man + delim_str \
  + figures_man + delim_str \
  + format_man \
  + begin_end_str
  
  # if no options are given
  if len(opts) == 1:
    sys.exit(man_string)
  
  for name in opts:
    if name.find('--help') >= 0 or name.find('-h') >= 0:
        sys.exit(man_string)

def check_for_figurestyles(opts):
  figstyles = ['b-', 'r-', 'g-']
  
  for name in opts:
    # if --format is given
    if name.find('--styles=') >= 0:
      cvs_styles = name.split('=')[1]
      figstyles = cvs_styles.split(',')
      break

  return figstyles

def check_for_figurenumbers(opts):
  fignumbers = [1, 3, 5, 6]
  
  for name in opts:
    # if --format is given
    if name.find('--figures=') >= 0:
      csv_figures = name.split('=')[1]
      figures = csv_figures.split(',')
      fignumbers = []
      for str in figures:
        fignumbers.append(int(str)) 
      break

  return fignumbers

def check_for_format(opts):
  fileformat = '*.out'
  
  for name in opts:
    # if --format is given
    if name.find('--format=') >= 0:
      fileformat = name.split('=')[1]
      break
   
  return fileformat

def check_for_endtime(opts):
  endtime = 15.
  
  for name in opts:
    # if --format is given
    if name.find('--endtime=') >= 0:
      endtime = float(name.split('=')[1])
      break
   
  return endtime

# main part
if __name__ == '__main__':

  # filenames to process
  # read from the arguments to the python script
  args = sys.argv
  # args[0] is always the calling script name
  
  check_for_help(args)
  
  fileformat = check_for_format(args)
  fignumbers = check_for_figurenumbers(args)
  
  # colors for the figure
  color = check_for_figurestyles(args)
  
  endtime=check_for_endtime(args)
     
  show_plots = True
  
  # was an argument passed to this script?
  # if not, complain
  if len(args) == 1:
    print 'no file or directory specified'
    sys.exit()
  else:
    filename = args[-1]
  
  
  if len(fignumbers) == 0:
    print 'No figures to plot specified.'
    sys.exit()

  if any( fn < 1 or fn > 6 for fn in fignumbers ):
    print 'Figure numbers must range between 1 and 6 inclusive.'
    sys.exit()
  
  if os.path.isdir(filename):
    args = [ args[0] ] + sorted(glob.glob(os.path.join(filename, fileformat)))
    print 'Found', str(len(args)-1), 'files matching format', fileformat
  elif not os.path.isfile(filename):
    if filename.find('=') >= 0:
      print 'Last argument must be a file or a directory.'
      sys.exit()
    else:
      print 'No such file or directory: \"', filename, '\"'
      sys.exit()

  ####
  # some general plotting parameters

  f_id = -1
  is_first_file = True
  is_last_file = False
  
  fig5_xaxis_label = 't'
  fig5_yaxis_label = r'$E_{kin}$'
  
  fig6_xaxis_label = fig5_xaxis_label
  fig6_yaxis_label = 'rel. vorticity thickness' 
  
  fontsize = 15
  fontsize1 = 12
  linewidth = 2
  markersize = 8
  
  #
  ###


  if len(color) < len(args)-1:
    print 'Not enough plotting styles provided.'
    more_cols = len(args)-len(color)
    cmap = plt.get_cmap('jet')
    color.append(cmap(np.linspace(0,1,more_cols)))

  for args_it in range(1, len(args)):
    f_id = f_id+1
    
    if(args_it == len(args)-1):
      is_last_file = True
    
    filename = args[args_it]
    # get the basename of the file
    base = os.path.basename(filename)
    
    print 'processing file:', base
  
    legend_str = base.split('_')
    legend_str = ' '.join(legend_str[5:7])
  
    starttime = dt = 0.
    timesteps_str = []
    timesteps = []
  
    with open(filename) as f:
      for line in f:
        # get the start time
        if line.startswith('STARTTIME:'):
          starttime = float(line.split()[1])
        #elif line.startswith('ENDTIME:'):
          #endtime = float(line.split()[1])
        elif line.startswith('TIMESTEPLENGTH:'):
          dt = float(line.split()[1])
        elif line.startswith('CURRENT TIME:'):
          timesteps_str.append( line.split()[-1] )
          timesteps.append( float( line.split()[-1] ) )
    
     
    numtimesteps = len(timesteps)
    
    # kinetic energy for every time step
    ekin = np.zeros(numtimesteps)
    
    # relative vorticity for every time step
    rel_vort = np.zeros(numtimesteps)
    
    # this one is fixed
    block_end_str = os.linesep
    
    ts_iter = 0
    ts_str = timesteps_str[ts_iter]
    
    block_begin_str = 'CURRENT TIME: ' + ts_str
    
    # how to look for kinetic energy
    ekin_str = ts_str + ' kinetic energy '
    
    # how to look for vorticity (only the second value is important)
    vort_str = ts_str + ' vorticity thickness'
    
    in_block = False
    
    data = []
    y_counter = 0
    x_counter = 0
    
    with open(filename) as f:
      for line in f:
        # if we are in a value block
        if in_block:
          # get ekin
          if line.startswith(ekin_str):
            ekin[ts_iter] = line.split()[-1]

          # get relative vorticity
          elif line.startswith(vort_str):
            rel_vort[ts_iter] = line.split()[-1]
            
          # check for end of block
          elif line.startswith(block_end_str):
            in_block = False
            
            # set new strings to look for (next time step)
            ts_iter = ts_iter + 1
            
            ts_str = timesteps_str[ts_iter]
    
            block_begin_str = 'CURRENT TIME: ' + ts_str
            ekin_str = ts_str + ' kinetic energy '
            vort_str = ts_str + ' vorticity thickness'
            
            # the actual data formatted as:
            # words[0] = 't'
            # words[1] = current time step
            # words[2] = 'y' or 'x'
            # words[3] = y or x coordinate
            # words[4] = 'mu' or 'mv'
            # words[5] = mean velocity in x or y direction
            # words[6] = 'rms_u' or 'rms_v'
            # words[7] = root mean square velocity in x or y direction
            # 
            # next two currently are ignored
            # words[8] = 'R11' or 'R22'
            # words[9] = values for R11 or R22
            # 
            # We just need this for the very last time step
          elif line.startswith('t '):
            if ts_iter == int((endtime-starttime)/dt)-1 :
              words = line.split()
              
              # label y or x on y axis of first plot
              if (words[2] == 'y'):
                fig1_xaxis_label = r'$\overline{u}$'
                fig1_yaxis_label = 'y'
                
                fig3_xaxis_label =  r'$\sqrt{ < u ^ 2 > }$'
                fig3_yaxis_label = fig1_yaxis_label
                
                y_counter = y_counter + 1
              else:
                fig2_xaxis_label = r'$\overline{v}$'
                fig2_yaxis_label = 'x'
                
                fig4_xaxis_label =  r'$\sqrt{ < v ^ 2 > }$'
                fig4_yaxis_label = fig2_yaxis_label
                
                x_counter = x_counter + 1
             
              
              # [ y or x, mu or mv, rms_u or rms_v ]              
              data.append( [ float(words[3]), 
                             float(words[5]),
                             float(words[7]) ] )
            
        # get the start time
        elif line.startswith(block_begin_str):         
          in_block = True
    
    data = np.array(data)
    
    '''
    PLOT
    '''
    
    # location
    location =  ['lower left', 'upper left', 'lower right', 'upper right']
    loc_no = 3;    
    
    '''
    PLOTTING PARAMETERS
    '''
    legend = legend_str
    ##############################################################
    #
    # Figure 1    y over mu
    #
    ##############################################################
    if any( fn == 1 for fn in fignumbers):
      if is_first_file:
        f1 = plt.figure(figsize=(8,6))
        ax1 = f1.add_subplot(111)
      
      datax = data[:y_counter,1]
      datay = data[:y_counter,0]
      
      ax1.plot(datax,datay,color[f_id],label = legend,
                   linewidth=linewidth,markersize=markersize,
                   fillstyle='none',markeredgewidth=1.5)
      
      # compute axis
      if is_first_file:
        f1minx = min(datax)
        f1maxx = max(datax)
        f1miny = min(datay)
        f1maxy = max(datay)
        
      else:
        # maybe axis changed, compute axis
        if f1minx  > min(datax):
            f1minx = min(datax)
        if f1maxx < max(datax):
            f1maxx = max(datax)
        if f1miny  > min(datay):
            f1miny = min(datay)
        if f1maxy < max(datay):
            f1maxy = max(datay)
      
      if is_last_file:          
        ax1.set_xlim(f1minx*1.1, f1maxx*1.1)
        ax1.set_ylim(f1miny, f1maxy)
                              
        ax1.tick_params(axis='both', labelsize=fontsize)
        ax1.set_xlabel(fig1_xaxis_label,fontsize=fontsize)
        ax1.set_ylabel(fig1_yaxis_label,fontsize=fontsize)
        ax1.legend(loc=location[loc_no],fontsize=fontsize1)
    
    ##############################################################
    #
    # Figure 2    x over mv
    #
    ##############################################################
    if any( fn == 2 for fn in fignumbers):
      if is_first_file:
        f2 = plt.figure(figsize=(8,6))
        ax2 = f2.add_subplot(111)
      
      datax = data[:y_counter:,1]
      datay = data[:y_counter:,0]
      
      ax2.plot(datax,datay,color[f_id],label = legend,
                   linewidth=linewidth,markersize=markersize,
                   fillstyle='none',markeredgewidth=1.5)
      
      # compute axis
      if is_first_file:
        f2minx = min(datax)
        f2maxx = max(datax)
        f2miny = min(datay)
        f2maxy = max(datay)
        
      else:
        # maybe axis changed, compute axis
        if f2minx  > min(datax):
            f2minx = min(datax)
        if f2maxx < max(datax):
            f2maxx = max(datax)
        if f2miny  > min(datay):
            f2miny = min(datay)
        if f2maxy < max(datay):
            f2maxy = max(datay)
      
      if is_last_file:       
        ax2.set_xlim(f2minx, f2maxx)
        ax2.set_ylim(f2miny, f2maxy)
      
        ax2.tick_params(axis='both', labelsize=fontsize)
        ax2.set_xlabel(fig2_xaxis_label,fontsize=fontsize)
        ax2.set_ylabel(fig2_yaxis_label,fontsize=fontsize)
        ax2.legend(loc=location[loc_no],fontsize=fontsize1)
      
    ##############################################################
    #
    # Figure 3    y over rms_u
    #
    ##############################################################
    if any( fn == 3 for fn in fignumbers):
      if is_first_file:
        f3 = plt.figure(figsize=(8,6))
        ax3 = f3.add_subplot(111)
      
      datax = data[:y_counter,2]
      datay = data[:y_counter,0]
      
      ax3.plot(datax,datay,color[f_id],label = legend,
                   linewidth=linewidth,markersize=markersize,
                   fillstyle='none',markeredgewidth=1.5)
      
      # compute axis
      if is_first_file:
        f3minx = min(datax)
        f3maxx = max(datax)
        f3miny = min(datay)
        f3maxy = max(datay)
        
      else:
        # maybe axis changed, compute axis
        if f3minx  > min(datax):
            f3minx = min(datax)
        if f3maxx < max(datax):
            f3maxx = max(datax)
        if f3miny  > min(datay):
            f3miny = min(datay)
        if f3maxy < max(datay):
            f3maxy = max(datay)
      
      if is_last_file:         
        ax3.set_xlim(f3minx*1.1, f3maxx*1.1)
        ax3.set_ylim(f3miny, f3maxy)
      
        ax3.tick_params(axis='both', labelsize=fontsize)
        ax3.set_xlabel(fig3_xaxis_label,fontsize=fontsize)
        ax3.set_ylabel(fig3_yaxis_label,fontsize=fontsize)
        ax3.legend(loc=location[loc_no],fontsize=fontsize1)
    
    ##############################################################
    #
    # Figure 4    x over rms_v
    #
    ##############################################################
    if any( fn == 4 for fn in fignumbers):
      if is_first_file:
        f4 = plt.figure(figsize=(8,6))
        ax4 = f4.add_subplot(111)
      
      datax = data[y_counter:,2]
      datay = data[y_counter:,0]
      
      ax4.plot(datax,datay,color[f_id],label = legend,
                   linewidth=linewidth,markersize=markersize,
                   fillstyle='none',markeredgewidth=1.5)
      
      # compute axis
      if is_first_file:
        f4minx = min(datax)
        f4maxx = max(datax)
        f4miny = min(datay)
        f4maxy = max(datay)
        
      else:
        # maybe axis changed, compute axis
        if f4minx  > min(datax):
            f4minx = min(datax)
        if f4maxx < max(datax):
            f4maxx = max(datax)
        if f4miny  > min(datay):
            f4miny = min(datay)
        if f4maxy < max(datay):
            f4maxy = max(datay)
      
      if is_last_file:       
        ax4.set_xlim(f4minx, f4maxx)
        ax4.set_ylim(f4miny, f4maxy)
      
        ax4.tick_params(axis='both', labelsize=fontsize)
        ax4.set_xlabel(fig4_xaxis_label,fontsize=fontsize)
        ax4.set_ylabel(fig4_yaxis_label,fontsize=fontsize)
        ax4.legend(loc=location[loc_no],fontsize=fontsize1)
    
    ##############################################################
    #
    # Figure 5    kinetic engery over time
    #
    ##############################################################
    if any( fn == 5 for fn in fignumbers):    
      if is_first_file:
        f5 = plt.figure(figsize=(8,6))
        ax5 = f5.add_subplot(111)
      
      datax = timesteps
      datay = ekin
      
      ax5.plot(datax,datay,color[f_id],label = legend,
                   linewidth=linewidth,markersize=markersize,
                   fillstyle='none',markeredgewidth=1.5)
      
      # compute axis
      if is_first_file:
        f5minx = min(datax)
        f5maxx = max(datax)
        f5miny = min(datay)
        f5maxy = max(datay)
        
      else:
        # maybe axis changed, compute axis
        if f5minx  > min(datax):
            f5minx = min(datax)
        if f5maxx < max(datax):
            f5maxx = max(datax)
        if f5miny  > min(datay):
            f5miny = min(datay)
        if f5maxy < max(datay):
            f5maxy = max(datay)
      
      if is_last_file:       
        ax5.set_xlim(f5minx, f5maxx)
        ax5.set_ylim(f5miny, f5maxy)
      
        ax5.tick_params(axis='both', labelsize=fontsize)
        ax5.set_xlabel(fig5_xaxis_label,fontsize=fontsize)
        ax5.set_ylabel(fig5_yaxis_label,fontsize=fontsize)
        ax5.legend(loc=location[loc_no],fontsize=fontsize1)
      
    ##############################################################
    #
    # Figure 6    rel vorticity thickness over time
    #
    ##############################################################
    if any( fn == 6 for fn in fignumbers):
      if is_first_file:
        f6 = plt.figure(figsize=(8,6))
        ax6 = f6.add_subplot(111)
      
      datax = timesteps
      datay = rel_vort
      
      ax6.plot(datax,datay,color[f_id],label = legend,
                   linewidth=linewidth,markersize=markersize,
                   fillstyle='none',markeredgewidth=1.5)
      
      # compute axis
      if is_first_file:
        f6minx = min(datax)
        f6maxx = max(datax)
        f6miny = min(datay)
        f6maxy = max(datay)
        
      else:
        # maybe axis changed, compute axis
        if f6minx  > min(datax):
            f6minx = min(datax)
        if f6maxx < max(datax):
            f6maxx = max(datax)
        if f6miny  > min(datay):
            f6miny = min(datay)
        if f6maxy < max(datay):
            f6maxy = max(datay)
      
      if is_last_file:
        ax6.set_xlim(f6minx, f6maxx)
        ax6.set_ylim(f6miny, f6maxy)
      
        ax6.tick_params(axis='both', labelsize=fontsize)
        ax6.set_xlabel(fig6_xaxis_label,fontsize=fontsize)
        ax6.set_ylabel(fig6_yaxis_label,fontsize=fontsize)
        ax6.legend(loc=location[loc_no],fontsize=fontsize1)
  
    is_first_file = False
  
  # end iteration over files
  #plt.tight_layout()
  plt.show()
