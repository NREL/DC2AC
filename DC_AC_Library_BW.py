# --------------------------------------------------------------------------------------------------
# DC-AC Tool functions
# Bin Wang
# 4/28/2022
#
# Reference:
# Wang, Bin, and Jin Tan. 2022. DC-AC Tool: Fully Automating the Acquisition of AC Power Flow
# Solution. Golden, CO: National Renewable Energy Laboratory. NREL/TP-6A40-80100.
# https://www.nrel.gov/docs/fy22osti/80100.pdf 
# --------------------------------------------------------------------------------------------------
from __future__ import with_statement
from contextlib import contextmanager


import glob, os, sys
import numpy
import csv
import xlrd



# --------------------------------------------------------------------------------------------------
# power flow data class
class PFData():
    def __init__(self):
        # bus data
        self.bus_num = []
        self.bus_type = []
        self.bus_Vpu = []

        # load data
        self.load_id = []
        self.load_bus = []
        self.load_Z = []
        self.load_I = []
        self.load_P = []
        self.load_MW = []
        self.load_Mvar = []
        self.load_MW_total = []
        self.load_Mvar_total = []

        # generator data
        self.gen_id = []
        self.gen_bus = []
        self.gen_status = []
        self.gen_S = []
        self.gen_mod = []
        self.gen_MW = []
        self.gen_Mvar = []
        self.gen_MW_total = []
        self.gen_Mvar_total = []

        # branch data
        self.brc_from = []
        self.brc_to = []
        self.brc_id = []
        self.brc_S = []
        self.brc_P = []
        self.brc_Q = []

    def getdata(self, psspy):
        # bus data
        self.bus_num = psspy.abusint(-1, 2, 'NUMBER')[1][0]
        self.bus_type = psspy.abusint(-1, 2, 'TYPE')[1][0]
        self.bus_Vpu = psspy.abusreal(-1, 2, 'PU')[1][0]

        # load data
        self.load_id = psspy.aloadchar(-1, 4, 'ID')[1][0]
        self.load_bus = psspy.aloadint(-1, 4, 'NUMBER')[1][0]
        self.load_Z = numpy.asarray(psspy.aloadcplx(-1, 4, 'YLACT')[1][0])
        self.load_I = numpy.asarray(psspy.aloadcplx(-1, 4, 'ILACT')[1][0])
        self.load_P = numpy.asarray(psspy.aloadcplx(-1, 4, 'MVAACT')[1][0])
        self.load_MW = self.load_Z.real + self.load_I.real + self.load_P.real
        self.load_Mvar = self.load_Z.imag + self.load_I.imag + self.load_P.imag
        self.load_MW_total = sum(self.load_MW)
        self.load_Mvar_total = sum(self.load_Mvar)

        # generator data
        self.gen_id = psspy.amachchar(-1, 4, 'ID')[1][0]
        self.gen_bus = psspy.amachint(-1, 4, 'NUMBER')[1][0]
        self.gen_status = psspy.amachint(-1, 4, 'STATUS')[1][0]
        self.gen_S = numpy.asarray(psspy.amachcplx(-1, 4, 'PQGEN')[1][0])
        self.gen_mod = numpy.asarray(psspy.amachint(-1, 4, 'WMOD')[1][0])
        self.gen_MW = self.gen_S.real
        self.gen_Mvar = self.gen_S.imag
        self.gen_MW_total = sum(self.gen_MW)
        self.gen_Mvar_total = sum(self.gen_Mvar)

        # branch data
        ierr, iarray = psspy.abrnint(-1, 0, 0, 3, 2, ['FROMNUMBER', 'TONUMBER'])
        self.brc_from = iarray[0][:]
        self.brc_to = iarray[1][:]
        self.brc_id = psspy.abrnchar(-1, 0, 0, 3, 2, ['ID'])[1][0]
        self.brc_S = numpy.asarray(psspy.abrncplx(-1, 1, 1, 3, 2, ['PQ'])[1][0])
        self.brc_P = self.brc_S.real
        self.brc_Q = self.brc_S.imag



# --------------------------------------------------------------------------------------------------
# check power flow solvability by running NR for up to multiple times (without adjusting loading)
# input: option - choose solver, n_flag1_max - maximum number of runs for power flow
# output: solved_flag, 0 - solved, 1 - reach maximum iteration number, 2 - blow up
def ifsolved(psspy, option, n_flag1_max):
    with silence():
        if option==1:
            psspy.fnsl([0, 0, 0, 1, 1, 0, 0, 0])
        else:
            psspy.fnsl([0, 0, 0, 1, 1, 1, 0])  # NR with flat start
            # psspy.fnsl([0, 0, 0, 1, 1, 1, 99, 0])  # first power flow with added generators & when checking solvability
        solved_flag = psspy.solved()

        n_flag1 = 0
        while solved_flag == 1:
            n_flag1 = n_flag1 + 1
            with silence():
                if option == 1:
                    psspy.nsol([0, 0, 0, 1, 1, 0, 0])   # fast decouple with "Do not flat start"
                else:
                    psspy.fnsl([0, 0, 0, 1, 1, 0, 0, 0])
                solved_flag = psspy.solved()
            if n_flag1 > n_flag1_max:
                print('    Flag = 1 does not disappear over ' + str(n_flag1_max) + ' runs. Let Flag = 2.')
                solved_flag = 2
    return solved_flag










def addgen(psspy, pfd, allsubs):
    # add generators to specified subs, specify P and V
    with silence():
        ad_pq_busid = []
        ad_pv_busid = []
        for ii in allsubs:
            bus_i = pfd.bus_num.index(ii)
            bus_i_type = pfd.bus_type[bus_i]
            # bus_i_Vpu = pfd.bus_Vpu[bus_i]

            if bus_i_type == 1:
                ad_pq_busid.append(ii)
                psspy.bus_chng_3(ii, [2, psspy._i, psspy._i, psspy._i],
                                 [psspy._f, 1.000, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f], psspy._s)
                psspy.plant_data(ii, psspy._i, [1.000, psspy._f])
                psspy.machine_data_2(ii, r"""ad""", [psspy._i, psspy._i, psspy._i, psspy._i, psspy._i, psspy._i],
                                     [psspy._f, psspy._f, 99999.0, -99999.0, 99999.0, -99999.0, psspy._f, psspy._f,
                                      psspy._f, psspy._f,
                                      psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f])
                # psspy.plant_chng(i, psspy._i, [bus_i_Vpu, 1])
            if bus_i_type == 2:
                ad_pv_busid.append(ii)
                psspy.bus_chng_3(ii, [2, psspy._i, psspy._i, psspy._i],
                                 [psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f], psspy._s)
                psspy.machine_data_2(ii, r"""ad""", [psspy._i, psspy._i, psspy._i, psspy._i, psspy._i, psspy._i],
                                     [psspy._f, psspy._f, 99999.0, -99999.0, 99999.0, -99999.0, psspy._f, psspy._f,
                                      psspy._f, psspy._f,
                                      psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f])
    return (ad_pq_busid, ad_pv_busid)



# solvability check by adjusting loading
def svbltcheck(psspy, pfd, step1_basecase, flag_raw_sav):
    # step 1: get a solvable case by reducing the loading
    dec_perc0 = float(0)
    dec_inc = float(5)
    dec_perc = dec_perc0
    n_iter = 0
    n_iter_max = 4
    solved_flag = 1

    insolvable_flag = 0
    while (solved_flag!= 0) & (n_iter < n_iter_max):
        n_iter = n_iter + 1
        dec_perc = dec_perc + dec_inc

        # early termination
        if dec_perc >= 10:
            insolvable_flag = 1
            break

        new_load_MW = (100 - dec_perc) / 100 * pfd.load_MW_total
        new_load_Mvar = (100 - dec_perc) / 100 * pfd.load_Mvar_total
        new_gen_MW = (100 - dec_perc) / 100 * pfd.gen_MW_total

        with silence():
            if flag_raw_sav == 1:
                psspy.read(0, step1_basecase)
            elif flag_raw_sav == 2:
                psspy.case(step1_basecase)

            psspy.scal_2(0, 1, 1, [0, 0, 0, 0, 0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
            psspy.scal_2(0, 1, 2, [psspy._i, 1, 0, 1, 0],
                         [new_load_MW, new_gen_MW, 0.0, -.0, 0.0, -.0, new_load_Mvar])

            psspy.fnsl([0, 0, 0, 1, 1, 1, 99, 0])
            solved_flag = psspy.solved()
    if insolvable_flag == 1:
        print('Power flow cannot converge above ' + str(
            100 - dec_perc) + '% loading.\n')
    return insolvable_flag, dec_perc


def approach_target_loading(psspy, pfd, dec_perc, step2_tempcase):
    inc_target = 100 / (100 - float(dec_perc))
    inc_cur = float(1.0001)
    inc_step = 0.01
    inc_step_MW = pfd.gen_MW_total / inc_target * inc_cur * inc_step
    inc_step_MW_min = float(1)

    pfdt = PFData()
    with silence():
        psspy.read(0, step2_tempcase)
    pfdt.getdata(psspy)

    for i in range(len(pfdt.gen_bus)):
        busi = pfdt.bus_num.index(pfdt.gen_bus[i])
        if pfdt.bus_type[busi] == 3 and pfdt.gen_status[i] == 1:
            pfdt.gen_MW_total = pfdt.gen_MW_total - pfdt.gen_MW[i]

    load_MW_pre = pfdt.load_MW_total
    gen_MW_pre = pfdt.gen_MW_total


    while (inc_cur < inc_target) & (inc_step_MW > inc_step_MW_min):
        inc_cur = inc_cur + inc_step
        with silence():
            psspy.read(0, step2_tempcase)


        # pfdt.getdata(psspy)

        gen_MW = gen_MW_pre * inc_cur
        load_MW = load_MW_pre * inc_cur
        # load_Mvar = load_Mvar_pre * inc_cur

        with silence():
            for i in range(len(pfdt.gen_bus)):
                busi = pfdt.bus_num.index(pfdt.gen_bus[i])
                if pfdt.bus_type[busi] == 3 and pfdt.gen_status[i] == 1:
                    psspy.machine_chng_2(pfdt.gen_bus[i],pfdt.gen_id[i],[psspy._i,psspy._i,psspy._i,psspy._i,psspy._i,psspy._i],[0.0,psspy._f,psspy._f,psspy._f,psspy._f,psspy._f,psspy._f,psspy._f,psspy._f,psspy._f,psspy._f,psspy._f,psspy._f,psspy._f,psspy._f,psspy._f,psspy._f])
            psspy.scal_2(0, 1, 1, [0, 0, 0, 0, 0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
            psspy.scal_2(0, 1, 2, [psspy._i, 1, 0, 1, 0], [load_MW, gen_MW, 0.0, -.0, 0.0, -.0, psspy._f])

        solved_flag = ifsolved(psspy, 1, 5)
        # pfdt.getdata(psspy)


        if solved_flag == 0:
            with silence():
                psspy.rawd_2(0, 1, [1, 1, 1, 0, 0, 0, 0], 0, step2_tempcase)
                # pfdt.getdata(psspy)
                # inc_cur_G = pfdt.gen_MW_total/load_MW_pre
                # inc_cur_G = inc_cur_G + inc_step
                pass
        else:
            inc_cur = inc_cur - inc_step
            inc_step = inc_step / 2
            # inc_cur_G = inc_cur_G - inc_step

        inc_step_MW = load_MW_pre * inc_step

        pass

    return inc_cur, inc_target, inc_step



def remove_added_gens(psspy, ad_pv_genbus, ad_pq_genbus, ad_pq_genQ_abs, step3_tempcase):
    # remove added gen on original pv bus
    n_ct_pv = 0
    for ii in ad_pv_genbus:
        with silence():
            psspy.read(0, step3_tempcase)
            psspy.purgmac(ii, r"""ad""")

        solved_flag = ifsolved(psspy, 1, 3)

        if solved_flag == 0:
            n_ct_pv = n_ct_pv + 1
            with silence():
                psspy.rawd_2(0, 1, [1, 1, 1, 0, 0, 0, 0], 0, step3_tempcase)
        if solved_flag == 2:
            break

    # remove added gen on original pq bus
    idx = sorted(range(len(ad_pq_genQ_abs)), key=lambda k: ad_pq_genQ_abs[k])
    n_ct_pq = 0
    n_step = 32
    n_cur = 0
    n_max = len(idx)
    while True:
        if n_step == 0:
            break
        if n_cur + n_step > n_max:
            n_step = n_max - n_cur

        with silence():
            psspy.read(0, step3_tempcase)
            for nn in range(n_cur, n_cur + n_step, 1):
                ii = idx[nn]
                psspy.purgmac(ad_pq_genbus[ii], r"""ad""")
                psspy.bus_chng_3(ad_pq_genbus[ii], [1, psspy._i, psspy._i, psspy._i],
                                 [psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f],
                                 psspy._s)

        solved_flag = ifsolved(psspy, 1, 3)

        if solved_flag == 0:
            n_ct_pq = n_ct_pq + n_step
            n_cur = n_cur + n_step
            with silence():
                psspy.rawd_2(0, 1, [1, 1, 1, 0, 0, 0, 0], 0, step3_tempcase)
        else:
            if n_step >= 1:
                if n_step == 2:
                    n_step = 1
                else:
                    n_step = int(n_step / 4)
            else:
                break
    return n_ct_pv, n_ct_pq, idx
    
    
    
    
    
    
    
def sol_out_orig(solved_flag):
    if solved_flag==0:
        print('Orig PF converged at target loading. (Success!)\n')
    else:
        print('Orig PF cannot converge at target loading.')

def read_subs(subs_filename):
    subs_file = xlrd.open_workbook(subs_filename)
    sheet = subs_file.sheet_by_index(0)
    N_subs = sheet.nrows

    allsubs = []
    for ii in range(N_subs):
        allsubs.append(int(sheet.cell_value(ii,0)))
    return allsubs


def removeCsvFile(outpath_csv, outpath_csv_Q, outpath_csv_V):
    if os.path.isfile(outpath_csv):
        os.remove(outpath_csv)
        # with open(outpath_csv, 'wb') as file:
        #     pass

    if os.path.isfile(outpath_csv_Q):
        os.remove(outpath_csv_Q)
        # with open(outpath_csv_Q, 'wb') as file:
        #     pass

    if os.path.isfile(outpath_csv_V):
        os.remove(outpath_csv_V)

def HaveOutputFolder(outpath, outpath_solved, outpath_solvedwQ, outpath_unsolved, outpath_temp, outpath_solved_Vadju):
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    if not os.path.exists(outpath_solved):
        os.mkdir(outpath_solved)
    if not os.path.exists(outpath_solvedwQ):
        os.mkdir(outpath_solvedwQ)
    if not os.path.exists(outpath_unsolved):
        os.mkdir(outpath_unsolved)
    if not os.path.exists(outpath_temp):
        os.mkdir(outpath_temp)
    if not os.path.exists(outpath_solved_Vadju):
        os.mkdir(outpath_solved_Vadju)





def saveRaw(psspy, rawfile, option, option2 = 0, perc = 0):
    outfile = str(rawfile)
    if option==1:
        outfile = outfile.replace("input", "dc2ac_output\solved")
        if option2 == 1:
            outfile = outfile[0:len(outfile) - 4] + """_solved.raw"""
        if option2 == 2:
            outfile = outfile[0:len(outfile) - 4] + """_addgen_remgen_solved.raw"""

    if option==2:
        outfile = outfile.replace("input", "dc2ac_output\\solvedwQ")
        outfile = outfile[0:len(outfile) - 4] + """_step3_addgen_solved.raw"""

    if option==3:
        outfile = outfile.replace("input", "dc2ac_output\unsolved")
        if option2 ==2:
            outfile = outfile[0:len(outfile) - 4] + """_step2_addgen_unsolved.raw"""
        if option2 ==2.1:
            outfile = outfile[0:len(outfile) - 4] + """_step2_addgen_solved_""" + perc + """.raw"""

    if option==4:
        outfile = outfile.replace("input", "temp")
        if option2==1:
            outfile = outfile[0:len(outfile) - 4] + """_step1_addgen.raw"""
        if option2==2:
            outfile = outfile[0:len(outfile) - 4] + """_step2_addgen_redP.raw"""
        if option2==2.1:
            outfile = outfile[0:len(outfile) - 4] + """_step2_addgen_solved.raw"""
        if option2==3:
            outfile = outfile[0:len(outfile) - 4] + """_step3_addgen_remgen.raw"""

    if option==5:
        outfile = outfile.replace("input", "dc2ac_output\\solved_Vadju")
        outfile = outfile[0:len(outfile) - 4] + """_step4_unlimQ.raw"""

    with silence():
        psspy.rawd_2(0, 1, [1, 1, 1, 0, 0, 0, 0], 0, outfile)

    return outfile


def getAddedGen(pfdt, ad_pq_busid, ad_pv_busid):
    ad_pq_genbus = []
    ad_pq_genQ = []
    ad_pq_genQ_abs = []

    ad_pv_genbus = []
    ad_pv_genQ = []

    for ii in range(len(pfdt.gen_bus)):
        if pfdt.gen_id[ii] == "AD":
            if pfdt.gen_bus[ii] in ad_pq_busid:
                ad_pq_genbus.append(pfdt.gen_bus[ii])
                ad_pq_genQ.append(pfdt.gen_Mvar[ii])
                ad_pq_genQ_abs.append(abs(pfdt.gen_Mvar[ii]))
            if pfdt.gen_bus[ii] in ad_pv_busid:
                ad_pv_genbus.append(pfdt.gen_bus[ii])
                ad_pv_genQ.append(pfdt.gen_Mvar[ii])
    return ad_pq_genbus, ad_pq_genQ, ad_pq_genQ_abs, ad_pv_genbus, ad_pv_genQ


def getQofAddedGen(psspy, step3_tempcase, ad_pq_busid, ad_pv_busid):
    with silence():
        psspy.read(0, step3_tempcase)
        psspy.fnsl([0, 0, 0, 1, 1, 0, -1, 0])
    temp_genid = psspy.amachchar(-1, 4, 'ID')[1][0]
    temp_genbus = psspy.amachint(-1, 4, 'NUMBER')[1][0]
    temp_genoutput = numpy.asarray(psspy.amachcplx(-1, 4, 'PQGEN')[1][0])
    temp_gen_MW = temp_genoutput.real
    temp_gen_Mvar = temp_genoutput.imag

    ad_genbus = []
    ad_genQ = []
    for ii in range(len(temp_genbus)):
        if temp_genid[ii] == "AD":
            if temp_genbus[ii] in ad_pq_busid:
                ad_genbus.append(temp_genbus[ii])
                ad_genQ.append(temp_gen_Mvar[ii])
            if temp_genbus[ii] in ad_pv_busid:
                ad_genbus.append(temp_genbus[ii])
                ad_genQ.append(temp_gen_Mvar[ii])
    return ad_genbus, ad_genQ


def writeResults(outpath_csv, outpath_csv_Q, outpath_csv_V, flag_summary, Q_support_summary, V_adjust_summary, rawfiles):
    with open(outpath_csv, 'wb') as file:
        writer = csv.writer(file)
        writer.writerow([" "])
        writer.writerow(
            ["Case names", "Solved directly", "Solved with Q support", "Solved with Vset adjusted", "Unsolved/Solved with reduced Pload"])
        n0 = 0
        n1 = 0
        n2 = 0
        n3 = 0
        for i in range(len(rawfiles)):
            rawfile = rawfiles[i]
            last_found = -1
            second_last = -1
            while True:
                second_last = last_found
                last_found = rawfile.find('\\', last_found + 1)
                if last_found == -1:
                    break
            idx = second_last
            if flag_summary[i] == 0:
                n0 = n0 + 1
                writer.writerow([rawfile[idx + 1:len(rawfile)], "1", "0", "0", "0"])
            if flag_summary[i] == 1:
                n1 = n1 + 1
                writer.writerow([rawfile[idx + 1:len(rawfile)], "0", "1", "0", "0"])
            if flag_summary[i] == 2:
                n2 = n2 + 1
                writer.writerow([rawfile[idx + 1:len(rawfile)], "0", "0", "1", "0"])
            if flag_summary[i] == 3:
                n3 = n3 + 1
                writer.writerow([rawfile[idx + 1:len(rawfile)], "0", "0", "0", "1"])
    print("Summary: (" + str(n0) + ", " + str(n1) + ", " + str(n2) + ", " + str(n3) + ")")

    with open(outpath_csv_Q, 'wb') as file:
        writer = csv.writer(file)
        writer.writerow([" "])
        writer.writerow(["Case names", "# of buses requiring Q support"])
        for i in range(len(rawfiles)):
            rawfile = rawfiles[i]
            last_found = -1
            second_last = -1
            while True:
                second_last = last_found
                last_found = rawfile.find('\\', last_found + 1)
                if last_found == -1:
                    break
            idx = second_last

            nn = int(Q_support_summary[i][0])
            stemp = str(Q_support_summary[i][0])
            for j in range(2 * nn):
                stemp = stemp + "," + str(Q_support_summary[i][j + 1])
            writer.writerow([rawfile[idx + 1:len(rawfile)], stemp])


    with open(outpath_csv_V, 'wb') as file:
        writer = csv.writer(file)
        writer.writerow([" "])
        writer.writerow(["Case names", "# of PV buses with adjusted Vset"])
        for i in range(len(rawfiles)):
            rawfile = rawfiles[i]
            last_found = -1
            second_last = -1
            while True:
                second_last = last_found
                last_found = rawfile.find('\\', last_found + 1)
                if last_found == -1:
                    break
            idx = second_last

            nn = int(V_adjust_summary[i][0])
            stemp = str(V_adjust_summary[i][0])
            for j in range(2 * nn):
                stemp = stemp + "," + str(V_adjust_summary[i][j + 1])
            writer.writerow([rawfile[idx + 1:len(rawfile)], stemp])




# By adjusting V set point at PV bus, those with close-to-zero Qinj can be removed directly (due to a better initial condition for NR)
def remove_geni_adjVi(psspy, rawfile, ad_genbus, ad_genQ, pfdt, step4_tempcase, delta_thrd = 1000):
    bus_remote = []
    Q_remote = []
    Vset_remote = []
    for busi, Qi in zip(ad_genbus, ad_genQ):
        Q = []
        with silence():
            psspy.read(0, step4_tempcase)
            psspy.plant_chng_4(busi, 0, [psspy._i, psspy._i], [1.00, psspy._f])
            flag = ifsolved(psspy, 1, 5)
            pfdt.getdata(psspy)
            idxi = pfdt.gen_bus.index(busi)
            while pfdt.gen_id[idxi] != 'AD':
                idxi = idxi + 1
            Q.append(pfdt.gen_Mvar[idxi])

            psspy.read(0, step4_tempcase)
            psspy.plant_chng_4(busi, 0, [psspy._i, psspy._i], [0.99, psspy._f])
            flag = ifsolved(psspy, 1, 5)
            pfdt.getdata(psspy)
            idxi = pfdt.gen_bus.index(busi)
            while pfdt.gen_id[idxi] != 'AD':
                idxi = idxi + 1
            Q.append(pfdt.gen_Mvar[idxi])

            psspy.read(0, step4_tempcase)
            psspy.plant_chng_4(busi, 0, [psspy._i, psspy._i], [1.01, psspy._f])
            flag = ifsolved(psspy, 1, 5)
            pfdt.getdata(psspy)
            idxi = pfdt.gen_bus.index(busi)
            while pfdt.gen_id[idxi] != 'AD':
                idxi = idxi + 1
            Q.append(pfdt.gen_Mvar[idxi])

            a = (-1.0000 * Q[0] + 0.5000 * Q[1] + 0.5000 * Q[2]) * 10000
            b = (2.0000 * Q[0] - 1.0050 * Q[1] - 0.9950 * Q[2]) * 10000
            c = (-0.9999 * Q[0] + 0.5050 * Q[1] + 0.4950 * Q[2]) * 10000

            delta = b * b - 4 * a * c

            if delta > delta_thrd:
                x1 = (-b + numpy.sqrt(b * b - 4 * a * c)) / 2 / a
                x2 = (-b - numpy.sqrt(b * b - 4 * a * c)) / 2 / a
                if abs(x1 - 1) < abs(x2 - 1):
                    vset = x1
                else:
                    vset = x2
                psspy.read(0, step4_tempcase)
                psspy.plant_chng_4(busi, 0, [psspy._i, psspy._i], [vset, psspy._f])
                psspy.fnsl([0, 0, 0, 1, 1, 0, 0, 0])
                flag = ifsolved(psspy, 1, 5)


                psspy.purgmac(busi, r"""ad""")
                psspy.bus_chng_3(busi, [1, psspy._i, psspy._i, psspy._i],
                                 [psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f],
                                 psspy._s)
                flag = ifsolved(psspy, 1, 3)
            else:
                flag = 2

            if flag == 0:
                step4_tempcase = saveRaw(psspy, rawfile, 5)
            else:
                print('cannot remove gen' + str(busi))
                psspy.read(0, step4_tempcase)

                pfdt.getdata(psspy)
                idxi = pfdt.gen_bus.index(busi)
                while pfdt.gen_id[idxi] != 'AD':
                    idxi = idxi + 1
                Qvset = pfdt.gen_Mvar[idxi]

                bus_remote.append(busi)
                Q_remote.append(Qvset)
                Vset_remote.append(vset)
    return Q_remote, bus_remote, Vset_remote


# get nearby buses - single layer
def getnearbybus(busi, pfdt):
    nearbus = []
    nearpv = []
    last_found = -1
    while busi in pfdt.brc_from[last_found + 1:]:
        last_found = pfdt.brc_from.index(busi, last_found + 1)
        if last_found == -1:
            break
        else:
            if (pfdt.brc_to[last_found] not in nearbus)&(pfdt.brc_to[last_found]!=busi):
                nearbus.append(pfdt.brc_to[last_found])
                idx = pfdt.bus_num.index(pfdt.brc_to[last_found])
                if pfdt.bus_type[idx] == 2:
                    nearpv.append(pfdt.brc_to[last_found])
    return nearbus, nearpv




# get nearby pv buses
def getnearpv(busi, pfdt, n_layer = 1):
    # get nearby buses and nearby pv buses
    n = 1
    nearbus, nearpv = getnearbybus(busi, pfdt)

    if len(nearpv) == 0:
        n_layer = n_layer + 1

    while n < n_layer:
        n = n + 1
        nearbus_pre = []
        for busii in nearbus:
            nearbus_pre.append(busii)

        for busii in nearbus_pre:
            nearbus_temp, nearpv_temp = getnearbybus(busii, pfdt)
            for ii in range(len(nearbus_temp)):
                if (nearbus_temp[ii] not in nearbus)&(nearbus_temp[ii]!=busi):
                    nearbus.append(nearbus_temp[ii])
            for ii in range(len(nearpv_temp)):
                if (nearpv_temp[ii] not in nearpv)&(nearpv_temp[ii]!=busi):
                    nearpv.append(nearpv_temp[ii])
        if len(nearpv) == 0:
            n_layer = n_layer + 1

    return nearpv, n_layer


def remove_geni_adjVj(busi, busj, psspy, pfdt, step5_tempcase, vseti):
    Q = []
    with silence():
        psspy.read(0, step5_tempcase)
        psspy.plant_chng_4(busi, 0, [psspy._i, psspy._i], [vseti, psspy._f])
        psspy.plant_chng_4(busj, 0, [psspy._i, psspy._i], [1.0, psspy._f])
        flag = ifsolved(psspy, 1, 5)
        pfdt.getdata(psspy)
        idxi = pfdt.gen_bus.index(busi)
        while pfdt.gen_id[idxi] != 'AD':
            idxi = idxi + 1
        Q.append(pfdt.gen_Mvar[idxi])

        psspy.read(0, step5_tempcase)
        psspy.plant_chng_4(busi, 0, [psspy._i, psspy._i], [vseti, psspy._f])
        psspy.plant_chng_4(busj, 0, [psspy._i, psspy._i], [1.01, psspy._f])
        flag = ifsolved(psspy, 1, 5)
        pfdt.getdata(psspy)
        idxi = pfdt.gen_bus.index(busi)
        while pfdt.gen_id[idxi] != 'AD':
            idxi = idxi + 1
        Q.append(pfdt.gen_Mvar[idxi])


        a = (-1.0000 * Q[0] + 1.0000 * Q[1]) * 100
        b = Q[0] - a
        if  a!= 0:
            vset = -b/a + 0.01

            psspy.read(0, step5_tempcase)
            psspy.plant_chng_4(busi, 0, [psspy._i, psspy._i], [vseti, psspy._f])
            psspy.plant_chng_4(busj, 0, [psspy._i, psspy._i], [vset, psspy._f])
            psspy.fnsl([0, 0, 0, 1, 1, 0, 0, 0])
            flag = ifsolved(psspy, 1, 5)

            # pfdt.getdata(psspy)

            psspy.purgmac(busi, r"""ad""")
            psspy.bus_chng_3(busi, [1, psspy._i, psspy._i, psspy._i],
                             [psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f],
                             psspy._s)
            flag = ifsolved(psspy, 1, 3)
        else:
            vset = 0
            flag = 2

    return vset, flag


# The following function was prepared based on this reference: Whit. 2012. “Silencing PSSE Output,” Python for Power Systems. [Online] http://www.whit.com.au/blog/2012/03/silencing-psse-output/ 
@contextmanager
def silence(file_object=None):
    """
    Discard stdout (i.e. write to null device) or
    optionally write to given file-like object.
    """
    if file_object is None:
        file_object = open(os.devnull, 'w')

    old_stdout = sys.stdout
    try:
        sys.stdout = file_object
        yield
    finally:
        sys.stdout = old_stdout