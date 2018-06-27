'''
a=0
for i in range(0,20):
    print "%d %d-%d"%(i,a,a+5)
    a+=5
a=0    
for i in range(20,30):
    print "%d %d-%d"%(i,a,a+10)
    a+=10
a=0
for i in range(30,35):
    print "%d %d-%d"%(i,a,a+20)
    a+=20
a=0
for i in range(35,37):
    print "%d %d-%d"%(i,a,a+40)
    a+=40
'''
#Cent_mat
a=0
b="string Cent_mat[]={"
c=[]
for i in range(0,20):
    b+='"%d-%d%%"'%(a,a+5)
    c.append('"%d-%d%%"'%(a,a+5))
    a+=5
    b+=","
b+="};"
print b
print "NCent_mat=%d"%len(c)
print "\n"

#Cent_mat_new
a=0
b="string Cent_mat_new[]={"
c=[]
for i in range(0,20):
    b+='"%d-%d%%"'%(a,a+5)
    c.append('"%d-%d%%"'%(a,a+5))
    a+=5
    b+=","
a=0
for i in range(20,30):
    b+='"%d-%d%%"'%(a,a+10)
    c.append('"%d-%d%%"'%(a,a+10))
    a+=10
    b+=","
a=0
for i in range(30,35):
    b+='"%d-%d%%"'%(a,a+20)
    c.append('"%d-%d%%"'%(a,a+20))
    a+=20
    b+=","
a=0
for i in range(35,37):
    b+='"%d-%d%%"'%(a,a+40)
    c.append('"%d-%d%%"'%(a,a+40))
    a+=40
    b+=","
b+="};"
print b
print "NCent_mat=%d"%len(c)
