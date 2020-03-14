from graphenemodeling.graphene import monolayer as mlg
import matplotlib.pyplot as plt


# Markers
mcolor='salmon'
mnpairs=[
	[0,0],
	[0,1],
	[1,1]
]

pos=[]
for i in range(len(mnpairs)):
	m = mnpairs[i][0]
	n = mnpairs[i][1]
	pos.append(mlg.AtomicPosition(m,n,0)*1e10)
	pos.append(mlg.AtomicPosition(m,n,1)*1e10)

pos.append(mlg.AtomicPosition(-1,1,1)*1e10)
pos.append(mlg.AtomicPosition(1,0,0)*1e10)

pos.append(mlg.AtomicPosition(2,1,0)*1e10)
pos.append(mlg.AtomicPosition(1,2,1)*1e10)


fig, ax = plt.subplots(figsize=(8,4))

def setfont(font):
    return r'\font\a %s at 14pt\a ' % font

y=-1
alpha=0.5
plt.text(0,y,setfont('sen')+'G',fontsize=250,alpha=alpha)
plt.text(2,y,'M',fontsize=250,alpha=alpha)

for i in range(len(pos)):
	ax.plot(pos[i][0],pos[i][1],".",
			markersize=20,color=mcolor)



ax.set_axis_off()
plt.show()