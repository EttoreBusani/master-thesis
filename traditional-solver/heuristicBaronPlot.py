import numpy as np
import matplotlib.pyplot as plt

nApps = 3
# Load data from output files
baronFile = "/mnt/d/EDGE_CLOUD_GAME/MEC_Baron/outputFiles/" + str(nApps) + "apps.txt"
#baronFile = "/mnt/c/Users/busan/OneDrive/Polimi/TESI/Ettore_Busani_10751402/outputFiles/heuristic" + str(4) + "apps.txt"
heuristicFile = "/mnt/c/Users/busan/OneDrive/Polimi/TESI/Ettore_Busani_10751402/outputFiles/heuristicNew" + str(nApps) + "apps.txt"

n_users, baronProfit, baronTime, baronEdgeNodes, baronCloudVMs, baronUsersEdge, baronUsersCloud = np.loadtxt(baronFile,unpack=True)
n_users, heuristicProfit, heuristicTime, heuristicEdgeNodes, heuristicCloudVMs, heuristicUsersEdge, heuristicUsersCloud, heuristicR = np.loadtxt(heuristicFile,unpack=True)

delta = (baronProfit-heuristicProfit)/baronProfit
for x in delta: print(x*100)
print('Final delta:')
print(sum(delta)/len(delta)*100)


# Set position of bar on X axis
barWidth = 0.3
br1 = np.arange(len(n_users))
br2 = [x + barWidth for x in br1]

cmap = plt.get_cmap("tab10")

#--------------------------------------------------------------------------------------------
# PROFIT & EXECUTION TIME
#--------------------------------------------------------------------------------------------

fontSize = 12
titleSize = 18
edgeNodes = n_users*0.3
edgeNodes = edgeNodes.astype(int)

fig3, ax3 = plt.subplots(1)
#fig2.suptitle('Performace Analysis: Heuristic vs Baron', fontsize = titleSize)

ax3.bar(br1, heuristicProfit, width = barWidth, label ='heuristic')
ax3.bar(br2, baronProfit, width = barWidth, label ='baron')
plt.sca(ax3)
plt.xlabel('n. users', fontsize = fontSize)
plt.ylabel('Platform profit', fontsize = fontSize)
plt.xticks([r + barWidth/2 for r in range(len(n_users))], n_users.astype(int))
plt.legend()
#plt.title('Platform Profit')
plt.legend()


fig4, ax4 = plt.subplots(1)
ax4.bar(br1, np.log10(heuristicTime), width = barWidth, label ='heuristic')
ax4.bar(br2, np.log10(baronTime), width = barWidth, label ='baron')
plt.sca(ax4)
plt.xlabel('n. users', fontsize = fontSize)
plt.ylabel('Execution time (sec) - log scale', fontsize = fontSize)
plt.xticks([r + barWidth/2 for r in range(len(n_users))], n_users.astype(int))
plt.legend()
#plt.title('Execution Time')
plt.legend()
plt.show()

#--------------------------------------------------------------------------------------------
# USED RESOURCES & USERS SERVED
#--------------------------------------------------------------------------------------------

fig1, ax1 = plt.subplots(1)

plt.sca(ax1)
ax1.bar(br1, heuristicEdgeNodes, width = barWidth, label ='heuristic edge nodes',
color = 'lightblue')
ax1.bar(br1, heuristicCloudVMs, bottom=heuristicEdgeNodes, width = barWidth,label = 'heuristic cloud VMs',
color = cmap(0))
ax1.bar(br2, baronEdgeNodes, width = barWidth, label ='baron edge nodes',
color = 'bisque')
ax1.bar(br2, baronCloudVMs, bottom=baronEdgeNodes, width = barWidth, label = 'baron cloud VMs',
color= cmap(1))
ax1.scatter((br1+br2)/2,edgeNodes, s = 20, color = 'gold', label = 'max edge nodes available', marker = 'o',zorder=20)
plt.xticks([r + barWidth/2 for r in range(len(n_users))], n_users.astype(int))
ax1.set_xlabel('n. users', fontsize = fontSize)
ax1.set_ylabel('Avg resources used', fontsize = fontSize)
ax1.legend()
#plt.title('Resource Allocation')

fig2, ax2 = plt.subplots(1)
ax2.bar(br1, heuristicUsersEdge, width = barWidth, label ='heuristic users edge',
color = 'lightblue')
ax2.bar(br1, heuristicUsersCloud, bottom=heuristicUsersEdge, width = barWidth,label = 'heuristic users cloud',
color = cmap(0))
ax2.bar(br2, baronUsersEdge, width = barWidth, label ='baron users edge',
color = 'bisque')
ax2.bar(br2, baronUsersCloud, bottom=baronUsersEdge, width = barWidth, label = 'baron users cloud',
color= cmap(1))
plt.sca(ax2)
ax2.set_xlabel('n. users', fontsize = fontSize)
ax2.set_ylabel('Avg number of users served', fontsize = fontSize)
ax2.legend()
plt.xticks([r + barWidth/2 for r in range(len(n_users))], n_users.astype(int))
#plt.title('Users allocation')

#fig1.suptitle('Users and Resource Allocation: Heuristic vs Baron', fontsize = titleSize)
plt.show()



# # Plot the data
# plt.plot(n_users, heuristicProfit)
# plt.show()
