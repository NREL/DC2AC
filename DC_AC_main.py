# --------------------------------------------------------------------------------------------------
# DC-AC Tool: main function
# Bin Wang
# 4/28/2022
#
# Reference:
# Wang, Bin, and Jin Tan. 2022. DC-AC Tool: Fully Automating the Acquisition of AC Power Flow
# Solution. Golden, CO: National Renewable Energy Laboratory. NREL/TP-6A40-80100.
# https://www.nrel.gov/docs/fy22osti/80100.pdf 
# --------------------------------------------------------------------------------------------------

import os, sys

sys_path_PSSE = r'C:\Program Files (x86)\PTI\PSSE34\PSSPY27'
sys.path.append(sys_path_PSSE)
os_path_PSSE = r'C:\Program Files (x86)\PTI\PSSE34\PSSBIN'
os.environ['PATH'] += ';' + os_path_PSSE
os.environ['PATH'] += ';' + sys_path_PSSE

os.environ['PATH'] += ';' + r'C:\Python27'

import psse34
import psspy
import natsort
from DC_AC_Library_BW import *




idef = psspy.getdefaultint()
rdef = psspy.getdefaultreal()



def main():
    # setting
    flag_raw_sav = 1  # 1 - power flow raw file, 2 - power flow sav file


    # path to power flow files
    if flag_raw_sav == 1:
        path = r"""casefile\input\*.raw"""
        rawfiles = [f for f in glob.glob(path)]
        pffiles = natsort.natsorted(rawfiles, reverse = False)
    elif flag_raw_sav == 2:
        path = r"""casefile\input\*.sav"""
        savfiles = [f for f in glob.glob(path)]
        pffiles = natsort.natsorted(savfiles, reverse=False)



    # path to pre-specified substation bus numbers
    subs_filename = r"""casefile\input\subs_bus.xlsx"""

    # paths to output files
    outpath = r"""casefile\dc2ac_output\\"""
    outpath_solved = r"""casefile\dc2ac_output\solved\\"""                #1
    outpath_solvedwQ = r"""casefile\dc2ac_output\solvedwQ\\"""            #2
    outpath_unsolved = r"""casefile\dc2ac_output\unsolved\\"""            #3
    outpath_temp = r"""casefile\temp\\"""                                 #4
    outpath_solved_Vadju = r"""casefile\dc2ac_output\solved_Vadju\\"""    #5
    outpath_csv = r"""casefile\dc2ac_output\summary.csv"""
    outpath_csv_Q = r"""casefile\dc2ac_output\summary_Q.csv"""
    outpath_csv_V = r"""casefile\dc2ac_output\summary_V.csv"""



    # remove existing .csv output files
    removeCsvFile(outpath_csv, outpath_csv_Q, outpath_csv_V)

    # read pre-specified substation bus numbers
    allsubs = read_subs(subs_filename)

    # make sure output folders exist
    HaveOutputFolder(outpath, outpath_solved, outpath_solvedwQ, outpath_unsolved, outpath_temp, outpath_solved_Vadju)



    # claim a few variables
    flag_summary = []   # summary of solvability
    pfd = PFData()      # original power flow data
    pfdt = PFData()     # temp power flow data
    Q_support_summary = numpy.zeros((len(pffiles), len(allsubs) * 2 + 1)) # numbers of solved with Q cases, and added buses
    V_adjust_summary = numpy.zeros((len(pffiles), len(allsubs) * 2 + 1))  # numbers of solved cases with Vset adjusted at some PV buses


    # process all cases
    for pffile, i in zip(pffiles, range(len(pffiles))):
        print(pffile)

        with silence():
            psspy.psseinit(50000)
            if flag_raw_sav == 1:
                psspy.read(0, pffile)
            elif flag_raw_sav == 2:
                psspy.case(pffile)

        pfd.getdata(psspy)  # get power flow data

        solved_flag = ifsolved(psspy, 1, 5)  ## Step 1
        sol_out_orig(solved_flag)


        # if directly solved, save the solved case and move to the next
        if solved_flag == 0:  ## IF-1, yes
            saveRaw(psspy, pffile, 1, 1)
            flag_summary.append(0)  ## Possibility a
            continue



        ## IF-1, no
        # try to solve the unsolved AC power flow
        # add generators
        ad_pq_busid, ad_pv_busid = addgen(psspy, pfd, allsubs)  ## Step 2

        # save 1st unsloved case
        saveRaw(psspy, pffile, 4, 1)


        # solve 1st power flow
        solved_flag_1st = ifsolved(psspy, 0, 5)  ## Step 3, IF-2

        # if 1st power flow not converging, check solvability
        if solved_flag_1st !=0:  ## IF-2, no
            print('First power flow (with added generators) cannot converge! Investigating solvability...')
            with silence():
                if flag_raw_sav == 1:
                    psspy.read(0, pffile)
                elif flag_raw_sav == 2:
                    psspy.case(pffile)
            insolvable_flag, dec_perc = svbltcheck(psspy, pfd, pffile, flag_raw_sav)  ## Step 4 - Step 6

            if insolvable_flag == 1: # still cannot solve  ## Possibility c,
                saveRaw(psspy, pffile, 3, 2)

                flag_summary.append(3)
                Q_support_summary[i][0] = -1
                continue

            # save the solved case
            step2_tempcase = saveRaw(psspy, pffile, 4, 2)  ## IF-3, yes

            # try to approach the target loading
            inc_cur, inc_target, inc_step = approach_target_loading(psspy, pfd, dec_perc, step2_tempcase)  ## Step 5


            # Recover the original loading and dispatch
            if inc_cur >= inc_target:
                solved_flag_1st = 0

                # recover loading
                with silence():
                    psspy.read(0, step2_tempcase)  # saveRaw outputs .raw files only
                    pfdt.getdata(psspy)

                    psspy.scal_2(0, 1, 1, [0, 0, 0, 0, 0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
                    psspy.scal_2(0, 1, 2, [psspy._i, 1, 0, 1, 0],
                                 [pfd.load_MW_total, pfdt.gen_MW_total/pfdt.load_MW_total*pfd.load_MW_total, 0.0, -.0, 0.0, -.0, psspy._f])
                    psspy.fnsl([0, 0, 0, 1, 1, 0, 0, 0])

                solved_flag = ifsolved(psspy, 1, 5)



                if solved_flag!=0:
                    print('Should have solved the case.\n')
                    continue
                else:
                    with silence():
                        psspy.rawd_2(0, 1, [1, 1, 1, 0, 0, 0, 0], 0, step2_tempcase)

                # recover dispatch
                if solved_flag==0:
                    for i in range(len(pfd.gen_MW)):
                        busi = pfd.bus_num.index(pfd.gen_bus[i])

                        if pfd.gen_status[i]!=1 or abs(pfd.bus_type[busi])!=2:
                            continue

                        with silence():
                            psspy.read(0, step2_tempcase)
                            psspy.machine_chng_2(pfd.gen_bus[i], pfd.gen_id[i], [psspy._i, psspy._i, psspy._i, psspy._i, psspy._i, psspy._i],
                                             [pfd.gen_MW[i], psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f])
                            psspy.fnsl([0, 0, 0, 1, 1, 0, 0, 0])
                            solved_flag = ifsolved(psspy, 1, 5)
                        if solved_flag == 0:
                            with silence.silence():
                                psspy.rawd_2(0, 1, [1, 1, 1, 0, 0, 0, 0], 0, step2_tempcase)  ## Possibility b
                        else:
                            print('Should not see this!!!\n')
                            pass



            else:
                inc_cur = inc_cur - inc_step
                perc = str(inc_cur / inc_target * 100)[0:7]
                with silence():
                    psspy.read(0, step2_tempcase)

                saveRaw(psspy, pffile, 3, 2.1, perc)

                print('PF (w added gens) not converge, but can converged at %' + perc + ' loading. (Failed!)\n')  ## Possibility c
                flag_summary.append(3)
                Q_support_summary[i][0] = -1

                continue



        ## IF-2, yes
        # if 1st power flow converges, try to remove added generators
        wQ_bus = []
        wQ_Q = []
        if solved_flag_1st==0:
            print('PF (w added gens) converged.')
            saveRaw(psspy, pffile, 4, 2.1)
            step3_tempcase = saveRaw(psspy, pffile, 4, 3)

            # try to gradually remove all added generators
            # gen data from PSSE
            pfdt.getdata(psspy)
            ad_pq_genbus, ad_pq_genQ, ad_pq_genQ_abs, ad_pv_genbus, ad_pv_genQ = getAddedGen(pfdt, ad_pq_busid, ad_pv_busid)


            print('Trying to remove added generators..')
            # remove added gen
            n_ct_pv, n_ct_pq, idx = remove_added_gens(psspy, ad_pv_genbus, ad_pq_genbus, ad_pq_genQ_abs, step3_tempcase)  ## Step 7


            # all added gens can be removed  ## IF-8
            if n_ct_pq == len(idx):  ## IF-8, yes
                print('Orig PF converged at target loading. (Success!)\n')
                saveRaw(psspy, pffile, 1, 2)  ## Possibility d
                flag_summary.append(0)
                continue



            ## IF-7, no
            # added gen cannot be removed completely
            # get Q at unremovable added gen
            ad_genbus, ad_genQ = getQofAddedGen(psspy, step3_tempcase, ad_pq_busid, ad_pv_busid)


            # voltage adjustment to remove more added gen (directly solved by NR with a better initial)
            step4_tempcase = saveRaw(psspy, pffile, 5)
            Q_remote, bus_remote, Vset_remote = remove_geni_adjVi(psspy, pffile, ad_genbus, ad_genQ, pfdt, step4_tempcase)


            if len(bus_remote)==0:
                print('Orig PF converged (w adjusted local Vset). (Success!)\n')
                saveRaw(psspy, pffile, 1, 2)  ## Possibility e
                flag_summary.append(2)
                continue

            # near-by PV bus voltage adjustment to remove even more added gen
            step5_tempcase = saveRaw(psspy, pffile, 5)
            pfdt.getdata(psspy)

            Vset = []
            Vbus = []
            for busi, Qi, vseti in zip(bus_remote, Q_remote, Vset_remote):
                nearpv, n_layer = getnearpv(busi, pfdt, 3)

                flag = 2
                for busj in nearpv:
                    GenOnline = 0
                    last_found = -1
                    while busj in pfdt.gen_bus[last_found + 1:]:
                        last_found = pfdt.gen_bus.index(busj, last_found + 1)
                        if last_found == -1:
                            break
                        else:
                            if pfdt.gen_status[last_found] == 1:
                                GenOnline = 1

                    if GenOnline == 1:
                        Vjset, flag_temp = remove_geni_adjVj(busi, busj, psspy, pfdt, step5_tempcase, vseti)
                    else:
                        continue

                    if flag_temp == 0:
                        flag = 0
                        Vbus.append(busj)
                        Vset.append(Vjset)
                        break

                if flag == 2:
                    wQ_bus.append(busi)
                    wQ_Q.append(Qi)




        if len(wQ_bus)==0:
            print('Orig PF converged (w adjusted remote Vset). (Success!)\n')
            saveRaw(psspy, pffile, 1, 2)  ## Possibility e

            V_adjust_summary[i][0] = len(Vbus)
            for ii in range(len(Vbus)):
                V_adjust_summary[i][2 * ii + 1] = Vbus[ii]
                V_adjust_summary[i][2 * ii + 2] = Vset[ii]

            flag_summary.append(2)
            continue

        else:
            print('PF converged only when w Q support. ' + str(
                len(wQ_bus)) + ' out of ' + str(
                len(ad_pq_genbus) + len(ad_pv_genbus)) + ' left. \n')
            saveRaw(psspy, pffile, 2)

            flag_summary.append(1)  ## Possibility f
            Q_support_summary[i][0] = len(wQ_bus)

            # output Q support location and amount
            # print(idx[n_ct_pq:])
            for ii in range(len(wQ_bus)):
                Q_support_summary[i][2 * ii + 1] = wQ_bus[ii]
                Q_support_summary[i][2 * ii + 2] = wQ_Q[ii]






    # summarize processing results and write results into csv files
    # print(rawfiles)
    # print(flag_summary)
    writeResults(outpath_csv, outpath_csv_Q, outpath_csv_V, flag_summary, Q_support_summary, V_adjust_summary, rawfiles)




# main function
main()








