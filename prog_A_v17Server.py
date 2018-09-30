import subprocess
import os
from easygui import *
import sys
from multiprocessing import Process
import re
import numpy
import random
import string
import math
from collections import OrderedDict
import shutil
from scipy.interpolate import *
import pdb
from threading import Thread
import time
import datetime
import socket

def fit_triples(trans_1,trans_2,trans_3,top_17,peaklist,file_num,A,B,C,DJ,DJK,DK,dJ,dK,sorted_full_list,thread_num):
    
    
    peak_list_1 = peaklist[0:int(len(peaklist)/4)]#splits peaks into 4 parts to speed up processing
    peak_list_2 = peaklist[int(len(peaklist)/4):int(len(peaklist)/2)]
    peak_list_3 = peaklist[int(len(peaklist)/2):int(len(peaklist)*0.75)]
    peak_list_4 = peaklist[int(len(peaklist)*0.75):len(peaklist)]

    p_1 = float(peak_list_1[0][0])
    p_2 = float(peak_list_1[-1][0])
    p_3 = float(peak_list_2[0][0])
    p_4 = float(peak_list_2[-1][0])
    p_5 = float(peak_list_3[0][0])
    p_6 = float(peak_list_3[-1][0])
    p_7 = float(peak_list_4[0][0])
    p_8 = float(peak_list_4[-1][0])
    
#    file_values=open("variabili.txt","w")
#    file_values.write(str(list_a))
#    file_values.write("\n\n")
#    file_values.write(str(list_b))
#    file_values.write("\n\n")
#    file_values.write(str(list_c))
#    file_values.write("\n\n")
#    file_values.write(str(peak_list_1) + "\n"+str(peak_list_2)+"\n"+str(peak_list_3)+"\n"+str(peak_list_4))
#    file_values.write("\n\n")
#    file_values.write(str(sorted_full_list))
#    file_values.write("\n\n")
#    file_values.write(str(top_17))
#    
#    file_values.close()
    
    final_omc = []
    triples_counter = 0
    output_file = ""
    regular_counter = 0
    #error_counter = 0

    for freq_1,inten_1,freq_2,inten_2,freq_3,inten_3,total_diff in sorted_full_list:
        peaks_triple= [(str(freq_1),str(inten_1)),(str(freq_2),str(inten_2)),(str(freq_3),str(inten_3))]
        
        input_file = ""
        input_file += "anisole                                         Wed Mar Thu Jun 03 17:45:45 2010\n"
        input_file += "   8  500   5    0    0.0000E+000    1.0000E+005    1.0000E+000 1.0000000000\n" # don't choose more than 497 check transitions or it will crash.
        input_file +="a   1  1  0  50  0  1  1  1  1  -1   0\n"
        input_file += "           10000  %s 1.0E+004 \n" % A
        input_file += "           20000  %s 1.0E+004 \n" % B
        input_file += "           30000  %s 1.0E+004 \n" % C
        input_file += "             200  %s 1.0E-025 \n" % DJ
        input_file += "            1100  %s 1.0E-025 \n" % DJK
        input_file += "            2000  %s 1.0E-025 \n" % DK
        input_file += "           40100  %s 1.0E-025 \n" % dJ
        input_file += "           41000  %s 1.0E-025 \n" % dK
        fh_par = open("default%s%s.par"%(str(file_num),str(thread_num)),'w')
        fh_par.write(input_file)
        fh_par.close()

        input_file = ""#the next part adds in the three peaks to be fit
        input_file += trans_1[2][0:2]+' '+trans_1[2][2:4]+' '+trans_1[2][4:6]+' '+\
                      trans_1[3][0:2]+' '+trans_1[3][2:4]+' '+trans_1[3][4:6]+'                      '+peaks_triple[0][0]+' 0.50 1.0000\n' 
        input_file += trans_2[2][0:2]+' '+trans_2[2][2:4]+' '+trans_2[2][4:6]+' '+\
                      trans_2[3][0:2]+' '+trans_2[3][2:4]+' '+trans_2[3][4:6]+'                      '+peaks_triple[1][0]+' 0.50 1.0000\n' 
        input_file += trans_3[2][0:2]+' '+trans_3[2][2:4]+' '+trans_3[2][4:6]+' '+\
                      trans_3[3][0:2]+' '+trans_3[3][2:4]+' '+trans_3[3][4:6]+'                      '+peaks_triple[2][0]+' 0.50 1.0000\n'
        counter = 0
        for line in top_17:#the hack that adds in the check transitions but doesn't use them in the fit
            input_file += line[2][0:2]+' '+line[2][2:4]+' '+line[2][4:6]+' '+\
                      line[3][0:2]+' '+line[3][2:4]+' '+line[3][4:6]+'                      '+'%s.0'%(str(counter))+' 0.00 1.0000\n'
            counter += 1
        fh_lin = open("default%s%s.lin"%(str(file_num),str(thread_num)), "w")
        fh_lin.write(input_file)
        fh_lin.close()        
        a = subprocess.Popen("./spfit%s default%s%s"%(str(file_num),str(file_num),str(thread_num)), stdout=subprocess.PIPE, shell=True)
        a.stdout.read()#used to let SPFIT finish

        const_list = []

        fh_var = open("default%s%s.var"%(str(file_num),str(thread_num)))
        for line in fh_var:
            if line[8:13] == "10000":
                temp_A = float(line[15:37])
                const_list.append("%.3f" %temp_A)
            if line[8:13] == "20000":
                temp_B = float(line[15:37])
                const_list.append("%.3f" %temp_B)
            if line[8:13] == "30000":
                temp_C = float(line[15:37])
                const_list.append("%.3f" %temp_C)
        
        fh_fit = open("default%s%s.fit"%(str(file_num),str(thread_num)))
        file_list = []
        for line in fh_fit:
            file_list.append(line)

        freq_list = []
        for x in range(len(file_list)):
            if file_list[-x][11:14] == "RMS":
                rms_fit = float(file_list[-x][22:32]) #note - assumes RMS fit error is less than 1 GHz.  Change 22 to 21 if this is a problem.
            if file_list[-x][5:6] == ":" and int(file_list[-x][3:5])>3:
                freq_list.append(file_list[-x][60:71])
            if file_list[-x][40:64]=="EXP.FREQ.  -  CALC.FREQ.":
                break
        read_fit = (const_list[0],const_list[1], const_list[2],freq_list)
        triples_counter +=1
        constants = read_fit[0:3]
        freq_17 = read_fit[3]
        freq_17.reverse()
        A_1 = float(constants[0])
        B_1 = float(constants[1])
        C_1 = float(constants[2])
        omc_list = []
        
        theor_inten_list = []
        for x in range(len(top_17)):
            temp_inten = 10**float(top_17[x][0])
            theor_inten_list.append(temp_inten)

        theor_inten_ratio = (max(theor_inten_list)/min(theor_inten_list))
        theor_inten_avg = (sum(theor_inten_list)/len(theor_inten_list))
        theor_inten_unitless_stdev = numpy.std(theor_inten_list)/theor_inten_avg

        
        for x in range(len(top_17)): #matches peaks in the top 17 to peaks in experimental peak list <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            qnum_up = top_17[x][2]
            qnum_low = top_17[x][3]
            real_omc = 1.0
            current_peak = float(freq_17[x])
            if current_peak>p_1 and current_peak<p_2:#conditionals to find proper portion of experimental peak list to loop through
                peaks_section = peak_list_1
                peak = p_2
                regular_counter +=1
            elif current_peak>p_3 and current_peak<p_4:
                peaks_section = peak_list_2
                peak = p_4
                regular_counter +=1
            elif current_peak>p_5 and current_peak<p_6:
                peaks_section = peak_list_3
                peak = p_6
                regular_counter +=1
            elif current_peak>p_7 and current_peak<p_8:
                peaks_section = peak_list_4
                peak = p_8
                regular_counter +=1
            elif current_peak>p_8 or current_peak<p_1:
                peaks_section = peaklist
                peak = p_8
                #error_counter +=1
                real_omc = 0.0# this is the omc if you throw out peaks that go over the edge of the spectrum
            else:
                peaks_section = peaklist
                peak = p_8
                regular_counter +=1
            old_omc = 100000.0
            for peak_freq,peak_inten in peaks_section: #find nearest peak in actual spectrum to the given top 20 peak
                omc = abs(current_peak-float(peak_freq))
                omc_low = abs(current_peak-float(peak))
                if omc>old_omc:
                    omc_low = old_omc
                    omc_inten = temp_inten                    
                    break
                old_omc = omc
                temp_inten = peak_inten
                omc_inten = temp_inten
            if real_omc == 1.0:
                real_omc = omc_low
                
            omc_list.append((omc_low, real_omc, omc_inten))# current_peak,qnum_up,qnum_low)) you can add in this extra output, but its slower
        omc_avg = [float(omc) for omc, real_omc, omc_inten in omc_list]
        real_omc_avg = [float(real_omc) for omc, real_omc, omc_inten in omc_list]
        omc_inten_scoring = [float(omc_inten) for omc, real_omc, omc_inten in omc_list]
        score = str(len([omc for omc in omc_avg if omc<2.0])) #scores the accuracy of the fit, currently based on a peak being within 2 MHz which may be too coarse
        avg = (sum(omc_avg)/len(omc_avg))+rms_fit
        real_avg = (sum(real_omc_avg)/len(real_omc_avg))+rms_fit

        omc_inten_ratio = (max(omc_inten_scoring)/min(omc_inten_scoring))
        omc_inten_avg = (sum(omc_inten_scoring)/len(omc_inten_scoring))
        omc_inten_unitless_stdev = numpy.std(omc_inten_scoring)/omc_inten_avg

        score_inten_penalty = 1

        if (omc_inten_ratio <= 0.5*theor_inten_ratio) or (omc_inten_ratio >= 1.5*theor_inten_ratio):
            score_inten_penalty += 1

        if (omc_inten_unitless_stdev <= 0.5*theor_inten_unitless_stdev) or (omc_inten_unitless_stdev >= 1.5*theor_inten_unitless_stdev):
            score_inten_penalty += 1
        
        penalized_avg = avg*score_inten_penalty

        if float(A_1)>=float(B_1) and float(B_1)>=float(C_1) and float(C_1)>0:
            if int(score)<10: #makes sorting work properly later
                score = '0'+score  
            output_file += 'score = '+' '+score+' '+"Const = "+str(A_1)+' '+str(B_1)+' '+str(C_1)+' '+"average omc = "+str(avg)+'  '+"avg w/out peaks over edge = "+str(real_avg)+' '+"avg w/ inten penalty = "+str(penalized_avg)+"\n"

            if real_avg <= 0.2: #appends good finds (RMS < 0.2 MHz, ignoring peaks over edge) to interim file for each processor
                interim_output = 'score = '+' '+score+' '+"Const = "+str(A_1)+' '+str(B_1)+' '+str(C_1)+' '+"average omc = "+str(avg)+'  '+"avg w/out peaks over edge = "+str(real_avg)+' '+"avg w/ inten penalty = "+str(penalized_avg)+"\n"
                fh_interim_good = open("interim_good_output%s.txt"%(str(file_num)), "a")
                fh_interim_good.write(interim_output)
                fh_interim_good.close()

            if triples_counter == 100000: #appends to file after every 100000 triples
                fh_final = open("final_output%s%s.txt"%(str(file_num),str(thread_num)), "a")
                fh_final.write(output_file)
                fh_final.close()
                triples_counter = 0
                output_file = ""
    
    fh_final = open("final_output%s%s.txt"%(str(file_num),str(thread_num)), "a")#writes separate file for each thread
    #print 'out of %s peaks there were %s peaks that werent in the experimental spectrum'%(regular_counter, error_counter) 
    fh_final.write(output_file)
    fh_final.close()
    


def thread_creation(list_a,list_b,list_c,trans_1,trans_2,trans_3,top_17,peaklist,file_num,A,B,C,DJ,DJK,DK,dJ,dK):
    
    
    temp_list_a = []
    temp_list_b = []
    temp_list_c = []

    pred_ratio = pow(10,max(float(trans_1[0]),float(trans_2[0]),float(trans_3[0]))-min(float(trans_1[0]),float(trans_2[0]),float(trans_3[0])))

    for freq_1,inten_1 in list_a:
      temp_diff = abs(float(trans_1[1])-float(freq_1))
      temp_list_a.append((freq_1,inten_1,temp_diff))

    for freq_2,inten_2 in list_b:
      temp_diff = abs(float(trans_2[1])-float(freq_2))
      temp_list_b.append((freq_2,inten_2,temp_diff))

    for freq_3,inten_3 in list_c:
      temp_diff = abs(float(trans_3[1])-float(freq_3))
      temp_list_c.append((freq_3,inten_3,temp_diff))

    unsorted_full_list = []

    for freq_1,inten_1,diff_1 in temp_list_a:
        for freq_2,inten_2,diff_2 in temp_list_b:
            for freq_3,inten_3,diff_3 in temp_list_c:
                real_ratio = max(float(inten_1),float(inten_2),float(inten_3))/min(float(inten_1),float(inten_2),float(inten_3))
                avg_diff = (diff_1+diff_2+diff_3)/3
                scaled_diff = avg_diff*(abs(real_ratio-pred_ratio)+1) # Freq. difference scaled by deviation of intensity ratio from predicted
                unsorted_full_list.append((freq_1,inten_1,freq_2,inten_2,freq_3,inten_3,scaled_diff))

    sorted_full_list = sorted(unsorted_full_list, key=lambda entry: entry[6])

    num_thread=4
    start=0
    stop = len(sorted_full_list)/num_thread
    t = []
    for num in range(num_thread):
        t.append(Thread(target=fit_triples, args=(trans_1,trans_2,trans_3,top_17,peaklist,file_num,A,B,C,DJ,DJK,DK,dJ,dK,sorted_full_list[start:stop],num)))
        start=stop
        stop=stop + len(sorted_full_list)/num_thread
        if num==(num_thread-2):
            stop=len(sorted_full_list)

   
    for thread in t:
        thread.start()
    for thread in t:
        thread.join()

    a = subprocess.Popen("cat final_output%s*.txt > final_output%s.txt"%(str(file_num),str(file_num)), shell=True)
    a.wait()
    
    os.system("sort -r 'final_output%s.txt'>sorted_final_out%s.txt"%(str(file_num),str(file_num)))#sorts output by score
    
if __name__ == '__main__':

    HOST = ''                 # Symbolic name meaning all available interfaces
    PORT = 50008              # Arbitrary non-privileged port
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.bind((HOST, PORT))
    s.listen(1)
    conn, addr = s.accept()
    print 'Connected by', addr
    
    input_parameter=str(sys.argv[1])
    

    trans_1 = []
    trans_2 = []
    trans_3 = []
    trans_1_peaks = []
    trans_2_peaks = []
    trans_3_peaks = []
    top_peaks_3cut = []
    peaklist = []
    
    file_param = open(input_parameter)
    
    for line in file_param:
        if line.split()[0] == "A:":
            A = line.split()[1] 
        if line.split()[0] == "B:":
            B = line.split()[1]
        if line.split()[0] == "C:":
            C = line.split()[1]
        if line.split()[0] == "DJ:":
            DJ = line.split()[1]
        if line.split()[0] == "DK:":
            DK = line.split()[1]
        if line.split()[0] == "DJK:":
            DJK = line.split()[1] 
        if line.split()[0] == "dJ:":
            dJ = line.split()[1]
        if line.split()[0] == "dK:":
            dK = line.split()[1]
        if line.split()[0] == "processors:":
            processors=int(line.split()[1])
        if line.split()[0] == "folder:":
            job_name = line.split()[1] + "S"
        if line.split()[0] == "num_server:":
            num_server = line.split()[1]
        if line.startswith("trans_1:"):
            for line1 in file_param:
                if line1==" \n":
                    break
                else:
                   trans_1.append(line1[0:len(line1)-1])
        if line.startswith("trans_2:"):
            for line2 in file_param:
                if line2==" \n":
                    break
                else:
                   trans_2.append(line2[0:len(line2)-1])
        if line.startswith("trans_3:"):
            for line3 in file_param:
                if line3==" \n":
                    break
                else:
                   trans_3.append(line3[0:len(line3)-1])
        if line.startswith("trans_1_peaks:"):
            for line4 in file_param:
                if line4==" \n":
                    break
                else:
                    value1, value2 = (line4[0:len(line4)-1]).strip('()').split(',')
                    trans_1_peaks.append((value1, value2))
        if line.startswith("trans_2_peaks:"):
            for line5 in file_param:
                if line5==" \n":
                    break
                else:
                    value1, value2 = (line5[0:len(line5)-1]).strip('()').split(',')
                    trans_2_peaks.append((value1, value2))
        if line.startswith("trans_3_peaks:"):
            for line6 in file_param:
                if line6==" \n":
                    break
                else:
                    value1, value2 = (line6[0:len(line6)-1]).strip('()').split(',')
                    trans_3_peaks.append((value1, value2)) 
        if line.startswith("top_peaks_3cut:"):
            for line7 in file_param:
                if line7==" \n":
                    break
                else:
                    top_peaks_3cut.append(line7[0:len(line7)-1])
        if line.startswith("peaklist:"):
            for line8 in file_param:
                if line8=="\n":
                    break
                else:
                    peaklist.append(line8[0:len(line8)-1])
    file_param.close()

        
    a = subprocess.Popen("mkdir %s"%(job_name), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    a.wait()
    
    y = subprocess.Popen("cp spfit %s/spfit"%(job_name), stdout=subprocess.PIPE, shell=True)
    y.stdout.read()   
    
    y = subprocess.Popen("cp spcat %s/spcat"%(job_name), stdout=subprocess.PIPE, shell=True)
    y.stdout.read() 
    
    os.chdir("%s"%job_name)
    
    f_time = open("time.txt","w")
    inizio = time.time()
    f_time.write("Start:" + str(datetime.datetime.now()) + "\n")
    
    trans_1=tuple(trans_1)
    trans_2=tuple(trans_2)
    trans_3=tuple(trans_3)
    top_peaks_clean = []
    for entry in top_peaks_3cut:
        clean = entry[2:43]
        re_split = clean.split("', '")
        tuples = tuple(re_split)
        top_peaks_clean.append(tuples)
    top_peaks_3cut = top_peaks_clean
    
    peaks_clean = []
    for entry in peaklist:
        re_split = entry.split(" ",1)
        tuples = list(re_split)
        peaks_clean.append(tuples)
    peaklist = peaks_clean
    
    
    
    for number in range(processors):
        y = subprocess.Popen("cp spfit spfit%s"%(number), stdout=subprocess.PIPE, shell=True)
        y.stdout.read()  
    new_list = []

    new_list = [(len(trans_1_peaks),"trans_1_peaks"),(len(trans_2_peaks),"trans_2_peaks"),(len(trans_3_peaks),"trans_3_peaks")]
    new_list.sort()
        
    list_key = []
    list_a_peaks = [vars()[new_list[0][1]],new_list[0][1]]
    list_b_peaks = [vars()[new_list[1][1]],new_list[1][1]]
    list_c_peaks = vars()[new_list[2][1]]

    random.shuffle(list_c_peaks)  # Shuffle so that each processor gets a range of values for the third peak, not processor 0 getting only the lowest frequencies.  

    list_c_list = []
    for num in range(processors):
        processors = float(processors)
        num = float(num)
        x = int((num)*(len(list_c_peaks)/processors))
        y = int(len(list_c_peaks)*((num+1)/processors))
        list_c_list.append(list_c_peaks[x:y])
    list_c_list.append("marker")
    vars()[new_list[0][1]] = list_a_peaks[0]
    vars()[new_list[1][1]] = list_b_peaks[0]
    vars()[new_list[2][1]] = list_c_list

    processors = int(processors)
    for num in range(processors):

        if trans_1_peaks[-1]=="marker":
            trans_x_peaks = trans_1_peaks[num]
            trans_y_peaks = trans_2_peaks
            trans_z_peaks = trans_3_peaks
        
        if trans_2_peaks[-1]=="marker":
            trans_x_peaks = trans_1_peaks
            trans_y_peaks = trans_2_peaks[num]
            trans_z_peaks = trans_3_peaks

        if trans_3_peaks[-1]=="marker":
            trans_x_peaks = trans_1_peaks
            trans_y_peaks = trans_2_peaks
            trans_z_peaks = trans_3_peaks[num]
        
        vars()["p%s"%str(num)] = Process(target=thread_creation, args=(trans_x_peaks,trans_y_peaks,trans_z_peaks,trans_1,trans_2,trans_3,top_peaks_3cut,peaklist,num,A,B,C,DJ,DJK,DK,dJ,dK))
    

    inizio_procc = time.time()
    f_time.write("Start procc:" + str(datetime.datetime.now()) + "\n")
    
    for num in range(processors):
        vars()["p%s"%str(num)].start()
    for num in range(processors):
        vars()["p%s"%str(num)].join()
    fine_procc=time.time()
    f_time.write("End procc:" + str(datetime.datetime.now()) + "\n")
    minutes, seconds = divmod(fine_procc-inizio_procc, 60)
    hours, minutes = divmod(minutes, 60)
    f_time.write("Time without socket:" + str(int(hours)) + ":" + str(int(minutes)) + ":" + str(int(seconds)) + "\n" ) 
    print "Fine procc: " + str(int(hours)) + ":" + str(int(minutes)) + ":" + str(int(seconds))
    
    a = subprocess.Popen('cat sorted_final_out*.txt |sort -t "=" -k 4 -n > sorted_omc_catS.txt', shell=True)
    a.wait()
    


    file_send = open("sorted_omc_catS.txt")
    data = file_send.read(1024)
    while (data):
        conn.send(data)
        data = file_send.read(1024)
#    else:
#        print "errore nell'invio della dimensione del file"
    
    fine=time.time()
    f_time.write("End:" + str(datetime.datetime.now()) + "\n")
    seconds_tot=fine-inizio
    minutes, seconds = divmod(seconds_tot, 60)
    hours, minutes = divmod(minutes, 60)
    f_time.write("Time tot:" + str(int(hours)) + ":" + str(int(minutes)) + ":" + str(int(seconds)) )
    print "Time tot proccess: "+str(int(hours)) + ":" + str(int(minutes)) + ":" + str(int(seconds)) 
    #print str(int(hours)) + ":" + str(int(minutes)) + ":" + str(int(seconds))

    conn.close()
    f_time.close()
    
    
    

    
    
    
    
    
    
    
    
    
