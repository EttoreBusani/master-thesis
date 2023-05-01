import numpy as np
import matplotlib.pyplot as plt

nApps = 2
nUsers = 20

sourceFile = "/mnt/c/Users/busan/OneDrive/Polimi/TESI/Ettore_Busani_10751402/outputFiles/"+str(nUsers)+"users"+str(nApps)+"appsProfitPlot.txt"

r, profit = np.loadtxt(sourceFile,unpack=True)
fontSize = 12

dropDepl = [6.67955, 6.29795, 7.87688, 8.66685, 9.14372]
dropProfit = [15.8062, 16.8272, 18.4356, 18.9303, 17.4866]

changeDepl = [6.38982, 7.60823]
changeProfit = [14.974, 18.1793]


# fig, ax = plt.subplots(1)
# plt.sca(ax)
# plt.scatter(r,profit, s=8)
# plt.scatter(dropDepl,dropProfit, color = 'gold', label = 'drop depl. points')
# plt.scatter(changeDepl,changeProfit, color = 'red', label = 'change depl. points')
# #plt.title('Platform profit vs price')
# ax.set_xlabel('Price fee (r)', fontsize = fontSize)
# ax.set_ylabel('Platform profit', fontsize = fontSize)
# plt.legend()
# plt.xticks([x for x in np.arange(6,10,1)])
#
# plt.show()
#
#
#
# ####################################################################################
# ####################################################################################
#
#
#
# # fileTabu = "/mnt/c/Users/busan/OneDrive/Polimi/TESI/Ettore_Busani_10751402/outputFiles/tabu.txt"
# #
# # p1,p2,p3,p4 = np.loadtxt(fileTabu,unpack=True)
# #
# # var = [10,20,30,40,50,60,70,80,90]
# #
# # fig, ax = plt.subplots(1, figsize = (10,6))
# # plt.sca(ax)
# # plt.plot(var,p1,'-o',label = 'Neighbor selection type 0, tabu list attribute type 0')
# # plt.plot(var,p2,'-x',label = 'Neighbor selection type 0, tabu list attribute type 1')
# # plt.plot(var,p3,'-s',label = 'Neighbor selection type 1, tabu list attribute type 0')
# # plt.plot(var,p4,'-^',label = 'Neighbor selection type 1, tabu list attribute type 1')
# # ax.set_xlabel('Variance of users', fontsize = fontSize)
# # ax.set_ylabel('Percentage increase in profit', fontsize = fontSize)
# # plt.xticks(var, ['10%','20%','30%','40%','50%','60%','70%','80%','90%'])
# # plt.legend()
# # plt.show()
#
# ####################################################################################
# ####################################################################################
#
folderPath = "/mnt/c/Users/busan/OneDrive/Polimi/TESI/Ettore_Busani_10751402/outputFiles/"
file1 = "timeHeuristic3apps.txt"
file2 = "timeHeuristic4apps.txt"
file3 = "timeHeuristic5apps.txt"
file4 = "timeHeuristic6apps.txt"
file5 = "timeHeuristic7apps.txt"

n,time1 = np.loadtxt(folderPath+file1,unpack=True)
n,time2 = np.loadtxt(folderPath+file2,unpack=True)
n,time3 = np.loadtxt(folderPath+file3,unpack=True)
n,time4 = np.loadtxt(folderPath+file4,unpack=True)
n,time5 = np.loadtxt(folderPath+file5,unpack=True)


fig, ax = plt.subplots(1)
plt.sca(ax)
plt.plot(n,np.log10(time2),'-x',label = '3 apps')
plt.plot(n,np.log10(time1),'-o',label = '4 apps')
plt.plot(n,np.log10(time3),'-s',label = '5 apps')
plt.plot(n,np.log10(time4),'-^',label = '6 apps')
plt.plot(n,np.log10(time5),'-*',label = '7 apps')

ax.set_xlabel('n. of users', fontsize = fontSize)
ax.set_ylabel('Time (sec)', fontsize = fontSize)
plt.xticks(np.arange(100,1100,100))
plt.legend()
plt.show()


####################################################################################
####################################################################################


# folderPath = "/mnt/d/EDGE_CLOUD_GAME/Ettore_Busani_10751402/outputFiles/"
#
# file1 = "profitHeuristic1perm.txt"
# file2 = "profitHeuristic2perm.txt"
# file3 = "profitHeuristic3perm.txt"
# file4 = "profitHeuristic4perm.txt"
#
#
# n,profit1 = np.loadtxt(folderPath+file1,unpack=True)
# n,profit2 = np.loadtxt(folderPath+file2,unpack=True)
# n,profit3 = np.loadtxt(folderPath+file3,unpack=True)
# n,profit4 = np.loadtxt(folderPath+file4,unpack=True)
#
#
# fig, ax = plt.subplots(1)
# plt.sca(ax)
#
# plt.plot(n,profit1,'-o',label = 'k=1')
# plt.plot(n,profit2,'-x',label = 'k=2')
# plt.plot(n,profit3,'-s',label = 'k=3')
# plt.plot(n,profit4,'-^',label = 'k=4')
#
#
# ax.set_xlabel('n. of users', fontsize = fontSize)
# ax.set_ylabel('Platform profit', fontsize = fontSize)
# plt.xticks(np.arange(100,1100,100))
# plt.legend()
# plt.show()
