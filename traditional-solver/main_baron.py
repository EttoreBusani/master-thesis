import pyomo.environ as pyo
from pyomo.opt import SolverFactory
import time
import multiprocessing as mpp
import os
from Model import PlatformModel

# Saves all results on instance
def solve_with_Baron(instance,k):
    # Initialize all variables
    instance.r = (instance.r_max+instance.r_min)/2
    for a in instance.Applications:
        instance.edge_VM_number[a]=0
        instance.cloud_VM_number[a]=0
        for i in instance.Users_per_Application[a]:
            instance.y_e[i]=0
            instance.y_c[i]=0
            instance.x_1[i]=0
            instance.x_2[i]=0
            instance.t_1[i]=0
            instance.t_2[i]=0
            instance.z[i]=0
    opt = SolverFactory("baron",executable='/home/busani/baron-lin64/baron')
    opt.options['threads'] = int(mpp.cpu_count())
    Path = "logFiles"
    if not os.path.exists(Path): os.makedirs(Path)
    start=time.time()
    nUsers = len(instance.Users)
    nApps = len(instance.Applications)
    logFileName = Path + "/" + str(nUsers)+"users"+str(nApps)+"apps"+str(k)
    opt.solve(instance,options={'MaxTime': 7200}, keepfiles=True,tee=True,logfile=logFileName)
    exec_time=time.time()-start
    profit = instance.Maximum_edge_profit("Value")
    # flag_exceeded = search_str(log_Baron, '*** Max. allowable time exceeded ***')
    # flag_infeasible = search_str(log_Baron, '*** Problem is infeasible ***')
    return profit, exec_time, instance

def printSol(instance):
    print("n. users: " + str(len(instance.Users)) + ", n. apps: " + str(len(instance.Applications)))
    print("Platform profit: " + str(instance.Maximum_edge_profit("Value")))
    print("r: " + str(instance.r("Value")))
    for a in instance.Applications:
        print(str(a) + ": n_e = " + str(instance.edge_VM_number[a]("Value")) + ", n_c = " + str(instance.cloud_VM_number[a]("Value")))
        for i in instance.Users_per_Application[a]:
            print("\t" + str(i) +": y_e = " + str(instance.y_e[i]("Value")) + ", y_c = " + str(instance.y_c[i]("Value")))

def write_profit(filename,nUsers,profit):
    with open(filename,'w') as f:
        f.write(str(nUsers) + " " + str(profit))


def main():
    nApps = 3
    outputFile = "outputFiles/" + str(nApps)+"apps.txt"
    df = open(outputFile,'w')

    folder = "configFiles/"
    n_instances = 5
    for n in range(10,40,5):
        finalProfit = 0
        finalTime = 0
        for k in range(0,n_instances-1):
            filename = folder + str(n)+"users" + str(nApps)+"apps" + str(k) + ".yaml"
            platform = PlatformModel()
            data, model = platform.create_platform_model(filename)
            instance = model.create_instance(data)
            # instance.pprint()
            profit, time, sol = solve_with_Baron(instance,k)
            finalProfit += profit/n_instances
            finalTime += time/n_instances
            printSol(sol)
            # Write on output file

        df.write(str(len(sol.Users)) + " " + str(finalProfit) + " " + str(finalTime))
        df.write("\n")

    return 0



if __name__ == "__main__":
    main()
