# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 15:47:27 2017
@author: Samorg
"""
#################################################################################################
#////////////////////////////////////////////Load Modules\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#################################################################################################
import time
tmps1=time.clock()
import sys
sys.path.append('G:\Education\Projects\SmartR\Python Scripts\SmartR\FinalScripts')
sys.path.append('G:\Education\Projects\SmartR\Python Scripts\SmartR')
from bson import ObjectId
from functions import *
import pymongo
import pprint
from detect_peaks import detect_peaks
from Pearsoncoef import *
from math import *
import pandas as pd
import textwrap
import csv
import numpy as np
from scipy.signal import savgol_filter
from math import *
import matplotlib.pyplot as plt
import cStringIO
##PDf Libraries 
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.lib.utils import ImageReader




def extentionAnalyse():
    global maxxx
    global minnn
    global exercise
    exercise="extention"         
    maxx=detect_peaks(filtreflex,mph=-50,mpd=20,edge='rising',valley='true')
    minn=detect_peaks(filtreflex,mph=40,mpd=30,edge='rising')
    k=0
    maxxx=[]
    minnn=[]
#################################################################################################
#////////////////////////////Filtrage minimal des mins et MAx\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#################################################################################################    
    for i in range(0,len(maxx-k)):
        if flex[maxx[i]]<70:
            maxxx.append(maxx[i])
    for i in range(0,len(minn)):
        if flex[minn[i]]>-20:
           minnn.append(minn[i])            
    minlen=len(minnn)
    mouvements=[]
    size=0
    maximus=np.zeros([len(maxxx),11])     
    maximus=getMaximus(maximus)      
    maximus=maximus[maximus[:, 6] <150]
    maximus=maximus[maximus[:, 6] >10]
    maximus=maximus[maximus[:, 3]-maximus[:, 0] > 15]
    maximus=maximus[maximus[:, 3]-maximus[:, 0] > 15]
    maximus=maximus[maximus[:, 3]-maximus[:, 0] < 150]
    notes=getNotes()
    Notes=pd.DataFrame(notes,columns=names)
           
           
           
    return

   
def squatAnalyse():
    global maxxx
    global minnn
    exercise="Squat"   
    maxx=detect_peaks(filtreflex,mph=30,mpd=40,edge='rising',show='true')
    minn=detect_peaks(filtreflex,mph=-60,mpd=20,edge='rising',valley='true',show='true')
    k=0
    maxxx=[]
    minnn=[]  
#################################################################################################
#////////////////////////////Filtrage minimal des mins et MAx\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#################################################################################################    
    for i in range(0,len(maxx-k)):
        if flex[maxx[i]]>35:
            maxxx.append(maxx[i])     
    for i in range(0,len(minn)):
        if flex[minn[i]]<30:
            minnn.append(minn[i]) 
    minlen=len(minnn)
    #associer chaque max a ses  deux minimums 
    mouvements=[]
    size=0
#################################################################################################
#////////////////////////////Separation des Mouvements    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#################################################################################################
    maximus=np.zeros([len(maxxx),11])              
    maximus=getMaximus(maximus)                         
    maximus=maximus[maximus[:, 0]-maximus[:, 3] > 15]
    maximus=maximus[maximus[:, 0]-maximus[:, 3] > 15]
    maximus=maximus[maximus[:, 0]-maximus[:, 3] < 150]
#################################################################################################
#//////////////////////////////    Notes des Mouvements   \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#################################################################################################

    return maximus


def Proprio_Analyse():
    exercise="prop"         
    global maximus
    gx=savgol_filter(hexdata[:,3],61,3); 
    gy=savgol_filter(hexdata[:,4],21,5);   
    gz=savgol_filter(hexdata[:,5],21,5);
    ff=savgol_filter(filtreflex,101,3)
    Min_Mvt=detect_peaks(ff,mph=-30,mpd=150,edge='rising',valley='true',show='true')
    Max_GyrX=detect_peaks(gx,mph=20,mpd=200,edge='rising',show='true')
    me=mean(gx[Max_GyrX])
    Max_GyrX=detect_peaks(gx,mph=me-20,mpd=70,edge='rising',show='true')
    Min_GyrX=detect_peaks(gx,mph=20,mpd=100,edge='rising',valley='true',show='true')
    me=mean(gx[Min_GyrX])
    Min_GyrX=detect_peaks(gx,mph=-me-20,mpd=100,edge='rising',valley='true',show='true')    
    maximus=maximus=np.zeros([len(Min_GyrX),11])
    #for i in range(0,len(Min_Mvt)-1):
        
    for i in range (len(Min_GyrX)):
        maximus[i][0]=searchmin_Avant(Min_GyrX[i],exercise,Min_Mvt)
        maximus[i][1]=Min_GyrX[i]
        maximus[i][2]=searchmin_Apres(Min_GyrX[i],Max_GyrX)      
        maximus[i][3]=searchmin_Apres(Min_GyrX[i],Min_Mvt) 
        maximus[i][4]=maximus[i][3]-maximus[i][0]
        maximus[i][5]=variationsignal(int(maximus[i][0]),int(maximus[i][3]),int(maximus[i][1]),'flex',filtreflex,'prop')          
    filtreMvtsProprio()    
 

    return maximus    
    
def getMaximus(maximus):
    global exercise
    global maxxx
    global minnn    
    for i in range (len(maximus)):        
        maximus[i][0]=flex[maxxx[i]]
        maximus[i][1]=maxxx[i]
        maximus[i][2]=searchmin_Avant(maxxx[i],exercise,minnn)
        maximus[i][3]=flex[int(maximus[i][2])]
        maximus[i][4]=searchmin_Apres(maxxx[i],minnn)
        maximus[i][5]=flex[int(maximus[i][4])]
        maximus[i][6]=maximus[i][4]-maximus[i][2]
    #if maximus[i][2] < maximus[i][1]:
    maximus=maximus[maximus[:, 6] != 0]
    maximus=maximus[maximus[:, 6] > 0]
#maximus=maximus[maximus[:, 6] < 200]
    maximus=maximus[maximus[:, 2] != 0]
    maximus=np.delete(maximus,trim_table(maximus),0)
   # maximus=maximus[maximus[:, 7] <200]    
    return maximus







##########################################################################################
#//////////////////////////////  PDf Generation    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
##########################################################################################

def exportPdf(Series):
#del c
    plt.ioff()
    fichier="C:/datatest/rapports/Rapport_%s"%x[0]["%s"%pp]["user"]["firstname"]+"_%s.pdf"%idd
    c = canvas.Canvas(fichier)
    c.drawImage("C:/datatest/panda15.jpg",500,750,70,70)
    c.setFontSize(20)
    
    c.drawString(210,800," Rapport De Séance ")
    c.setFontSize(13)
    c.drawString(20,750," Nom:%s"%x[0]["%s"%pp]["user"]["firstname"])
    c.drawString(20,730," Prenom:%s"%x[0]["%s"%pp]["user"]["firstname"])
    c.drawString(20,710," Date:%s"%str(x[0]["%s"%pp]["createdAt"].ctime()))
    c.drawString(20,690," physiotherapist:%s"%x[0]["%s"%pp]["user"]["physiotherapist"])
    c.drawString(20,670," Exercice:%s"%x[0]["%s"%pp]["exercise"])
    
    if "comment" in x[0]["%s"%pp]:
        c.drawString(20,650," Commentaire:%s"%x[0]["%s"%pp]["comment"])
    else :
        c.drawString(20,650," Pas de commentaires ")
        
    
    #for i in range(len(Series)):
    #    c.drawString(100,700-i*15,"Bruit pour chaque mouvement :(%d)"%i+"est de(%d)" %maximus[i][7])
    for i in range(0,len(Series)):
        c.addPageLabel(i)
        fig = plt.figure(figsize=(20, 10))
        imgdata = cStringIO.StringIO()
        fig.suptitle('Angle de flexion', fontsize=20)
        plt.xlabel('temps ', fontsize=18)
        plt.ylabel('amplitude', fontsize=16)
        plt.plot(filtreflex[int(Series[i][0]):int(Series[i][1])])
        fig.savefig(imgdata, format='png')
        imgdata.seek(0)
        Image = ImageReader(imgdata)
       # c.drawString(250,650,"Angle de Flexion")
        #c.drawAlignedString(250,670,"Angle de Flexion")
        if i==0:
            c.drawImage(Image, 50, 390, 500, 250,preserveAspectRatio=True)
        else:
            c.drawImage(Image, 50, 520, 500, 250,preserveAspectRatio=True)
        
        fig2 = plt.figure(figsize=(20, 10))
        imgdata2= cStringIO.StringIO()
        plt.plot(filtrerot[int(Series[i][0]):int(Series[i][1])],color='green')
        fig2.suptitle('Angle de Rotation ', fontsize=20)
        plt.xlabel('temps ', fontsize=18)
        plt.ylabel('amplitude', fontsize=16)
        fig2.savefig(imgdata2, format='png')
        imgdata2.seek(0)
        Image2 = ImageReader(imgdata2)
        #c.drawString(250,380,"Angle de Rotation")
        if i==0:
            c.drawImage(Image2, 50, 100, 500, 250,preserveAspectRatio=True)  
        else:
            c.drawImage(Image2, 50, 170, 500, 250,preserveAspectRatio=True) 
        c.save()
        del(fig2)
##########################################################################################
#//////////////////////////////  PDf Generation    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
##########################################################################################
    


def getNotes():
    notes=np.zeros([len(maximus),11])
    for i in range(len(maximus)):
        maximus[i][7]=variationsignal(int(maximus[i][2]),int(maximus[i][4]),int(maximus[i][1]),'flex',filtreflex,exercise)          
        maximus[i][8]=Symetrieflex(int(maximus[i][2]),int(maximus[i][1]),int(maximus[i][4]),filtreflex)
        maximus[i][9]=np.float64(SymetrieRot(int(maximus[i][2]),int(maximus[i][1]),int(maximus[i][4]),filtrerot))
        maximus[i][10]=variationsignal(int(maximus[i][2]),int(maximus[i][4]),int(maximus[i][1]),'rot',filtrerot,exercise)      
        notes[i][0]=maximus[i][1]
        #note Bruit
        notes[i][1]=1-(maximus[i][7]/maximus[i][6])
        notes[i][2]=1-(maximus[i][10]/maximus[i][6])
        #note Symetrie Flex
        notes[i][3]=abs(maximus[i][8])
        #note symetrie Rot
        notes[i][4]=abs(maximus[i][9])   
        #correlation
        notes[i][5]=Corr_flex_rot(int(maximus[i][2]),int(maximus[i][1]),int(maximus[i][4]),exercise,filtrerot,filtreflex)/(maximus[i][4]-maximus[i][2])
        #note Max
        notes[i][6]=noterMax(notes[i][0],Ref,exercise,filtreflex)
        #bruit fourier
        #print(i)
        notes[i][7]=fourierTransform(int(maximus[i][1]),int(maximus[i][4]),'flex',notes[i][0],filtreflex,filtrerot)
        notes[i][8]=fourierTransform(int(maximus[i][1]),int(maximus[i][4]),'rot',notes[i][0],filtrerot,filtrerot)    
        notes[i][9]=(notes[i][1]*2+notes[i][3]*2+notes[i][6]*3+notes[i][7])/8
        notes[i][10]=(notes[i][2]+notes[i][4]+notes[i][8])/3
    return notes








def filtreMvtsProprio():
    global maximus
    todelete=[]
    k=0
    bo=False
    bo1=True
    i=0
    kk=0
    while bo1:
        kk+=1
        
        
        if i<=len(maximus)-2:    
            #print(i)
            if maximus[i][0]==maximus[i+1][0] :
                for  j in range(i,len(maximus)-1):
                    print(j)
                    if maximus[j][0]==maximus[j+1][0] :
                        todelete.append(j)
                        i=j
                        k=j
                    elif j==len(maximus)-2:
                       # print("sssss")
                        bo=True
                        break
            else:
                i=i+1                
                    
        if bo or kk>50:
            bo1=False
            break  
    
    maximus=np.delete(maximus,todelete,0)
    maximus=maximus[maximus[:, 4] != 0]
    maximus=maximus[maximus[:, 4] > 0]    
    
    
    
    
    
    
#################################################################################################
#//////////////////////////////   Detection Des Series    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#################################################################################################

def detectSeries():
    global maximus
    meanfrq=[]
    for i in range(1,len(maximus)-1):
        meanfrq.append(maximus[i][1]-maximus[i-1][1])
    
    
    Series=np.zeros([len(maximus),2])
    l=0    
    i=0
    b=True
    s=0
    m=0
    #minfrq=180
    minfrq=mean(meanfrq)
    while i <len(maximus) and b==True:
        s=s+1
        m=i+1
        #if s>len(maximus):
            #b=False
        dist=0
        dist2=0
        k=0
        for j in range(m,len(maximus)):
            #↓print(j)
            dist2=maximus[j-1][1]-maximus[i][1]
            #Sdist3=maximus[j][1]-maximus[i][1]
            dist=maximus[j][1]-maximus[i][1]
           # print("i",i,"j",j)
            #print("d1",dist,"d2",dist2)
            
            #print(dist-dist2)
            if (dist-dist2)<minfrq+30:
                k+=1
                #print(k)
            else:
                if k>=4:
                    Series[l][0]=maximus[i][2]
                    Series[l][1]=maximus[i+4][4]  
                    l+=1
                    i=j       
                    break
                elif k<4:
                  
                    Series[l][0]=maximus[i][2]  
                    Series[l][1]=maximus[i+k][4]
                    l+=1
                   
                    i=j
                    break
            
            if j==(len(maximus)-1) and i!=0 and k!=0:
                if k>=4:
                    Series[l][0]=maximus[i][2]  
                    Series[l][1]=maximus[i+4][4]
                    l+=1
                elif k<4:
                    #print(i)
                    #print(l)               
                    Series[l][0]=maximus[i][2]  
                    Series[l][1]=maximus[i+k][4]
                    l+=1
                    break
            
            if j==(len(maximus)-1) and i==0 and k!=0:
                if len(maximus)>4 and k>=4:
                    Series[0][0]=maximus[i][2]  
                    Series[0][1]=maximus[i+4][4]
                    l+=1
                elif len(maximus)>4 and k<4:
                    #∟print(i)
                    #print(l)               
                    Series[0][0]=maximus[i][2]  
                    Series[0][1]=maximus[i+k][4]
                    l+=1
                elif len(maximus)<4 :                
                    Series[0][0]=maximus[i][2]  
                    Series[0][1]=maximus[i+len(maximus)-1][4]
                    l+=1                
                                
                b=False        
                break                  
        if s>80:
            b=False
    return Series            
            
def detectSeriesProp():
    meanfrq=[]    
    for i in range(1,len(maximus)-1):
        meanfrq.append(maximus[i][0]-maximus[i-1][0])
            
    Series=np.zeros([len(maximus),2])
    l=0    
    i=0
    b=True
    s=0
    m=0
    #minfrq=100
    minfrq=mean(meanfrq)
    while i <len(maximus) and b==True:
        s=s+1
        m=i+1
        #if s>len(maximus):
            #b=False
        dist=0
        dist2=0
        k=0
        for j in range(m,len(maximus)):
            #↓print(j)
            dist2=maximus[j-1][0]-maximus[i][0]
            #Sdist3=maximus[j][1]-maximus[i][1]
            dist=maximus[j][0]-maximus[i][0]
           # print("i",i,"j",j)
            #print("d1",dist,"d2",dist2)
            
            #print(dist-dist2)
            if (dist-dist2)<minfrq+minfrq/2:
               # print(k)
                
                k+=1
            else:
                if k>=4:
    
                    Series[l][0]=maximus[i][0]
                    Series[l][1]=maximus[i+4][3]  
                    l+=1
                    i=j       
                    break
                elif k<4:
                  
                    Series[l][0]=maximus[i][0]  
                    Series[l][1]=maximus[i+k][3]
                    l+=1
                   
                    i=j
                    break
            
            if j==(len(maximus)-2) and i!=0 and k!=0:
                print(j)
                if k>=4:
                    Series[l][0]=maximus[i][0]  
                    Series[l][1]=maximus[i+4][3]
                    l+=1
                elif k<4:
                    #∟print(i)
                    print(k)               
                    Series[l][0]=maximus[i][0]  
                    Series[l][1]=maximus[i+k][3]
                    l+=1
                    break
                
            if j==(len(maximus)-1) and i==0 and k!=0:
                if len(maximus)>4 and k>=4:
                    Series[0][0]=maximus[i][0]  
                    Series[0][1]=maximus[i+4][3]
                    l+=1
                elif len(maximus)>4 and k<4:
                    #∟print(i)
                    #print(l)               
                    Series[0][0]=maximus[i][0]  
                    Series[0][1]=maximus[i+k][3]
                    l+=1
                elif len(maximus)<4 :                
                    Series[0][0]=maximus[i][0]  
                    Series[0][1]=maximus[i+len(maximus)-1][3]
                    l+=1                
                    
                    break
            
                b=False        
                break 
        if s>80:
            b=False
    
    


    return Series
    
#################################################################################################
#////////////////////////////////////////////Db_Data\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#################################################################################################









idd=ObjectId("5981dde3ff77e008ed614f0d")
MONGO_HOST = "127.0.0.1"
MONGO_DB = "Panda2"
names=['id','bruitFlex','bruitRot','symflex','symrot','cor','Note MAx','fftFlex','fftRot','Noteflex','NoteRot']
Ref=110
client = pymongo.MongoClient(MONGO_HOST,27017)
db = client[MONGO_DB]
exer="squats_1foot"

collection = db.exer2
#result = collection.find({'$and':[{"exercise":exer},{'user.firstname':fname},{"_id":idd}]})
result = collection.find({"_id":idd})

x = []
for i in result:
    x.append(i)
#print(x)
d=[]
spliteddata=[]
decrypteddata=[]
pp=9

#len(x[0])
Data=x[0]["%s"%pp]["messages"]
idd=x[0]["%s"%pp]["_id"]
exer=x[0]["%s"%pp]["exercise"]
fields = Data.split('],')
split1,split4,split3,split2=[],[],[],[]
for i in range(len(fields)):
            #print(i)
    if i==0:
        val = fields[i].split('[[', 1)[1].split(']')[0]
        #print(val)
    else:
        val = fields[i].split('[', 1)[1].split(']')[0]
               # print(val)
    d.append(val)
    spliteddata.append(d[i].split(','))
hexdata=np.zeros([len(spliteddata),17])
spliteddata2=[]
for i in range(0,len(spliteddata)):
    spliteddata2.append(spliteddata[i][1:21])
#Shift and split  of bytes    
    split1.append(splitbits(spliteddata2[i][4]))
    split2.append(splitbits(spliteddata2[i][9]))
    split3.append(splitbits(spliteddata2[i][14]))
    split4.append(splitbits(spliteddata2[i][19]))
#Bitwise operations     
    for m in range(0,4):
        hexdata[i][m]=int(spliteddata2[i][m])|split1[i][m]<<8 
    for m in range(4,8):
        hexdata[i][m]=int(spliteddata2[i][m+1])|split2[i][m-4]<<8       
    for m in range(8,12):
        hexdata[i][m]=int(spliteddata2[i][m+2])|split3[i][m-8]<<8
    for m in range(12,16):
        hexdata[i][m]=int(spliteddata2[i][m+3])|split4[i][m-12]<<8      
               # d2.append(d[0][i].split(",",20))
    hexdata[i][0]=hexdata[i][0]/10 -50    
    hexdata[i][1]=hexdata[i][1]/10 -50    
    hexdata[i][2]=hexdata[i][2]/10 -50    
    hexdata[i][3]=hexdata[i][3]-500    
    hexdata[i][4]=hexdata[i][4]-500    
    hexdata[i][5]=hexdata[i][5]-500    
    hexdata[i][6]=hexdata[i][6]/10 -50    
    hexdata[i][7]=hexdata[i][7]/10 -50    
    hexdata[i][8]=hexdata[i][8]/10 -50    
    hexdata[i][9]=hexdata[i][9]-500    
    hexdata[i][10]=hexdata[i][10]-500    
    hexdata[i][11]=hexdata[i][11]-500    
    hexdata[i][12]=(hexdata[i][12]-200)/10     
    hexdata[i][13]=(hexdata[i][13]-200)/10    
    hexdata[i][14]=hexdata[i][14]
    hexdata[i][15]=hexdata[i][15]   
    hexdata[i][16]=float(spliteddata[i][0])   

    i=i+1          
    
data=hexdata
del split1,split2,split3,split4
flexion=np.zeros([len(data),3])
rotation=np.zeros([len(data),4])










#################################################################################################
#////////////////////////////Angles Calculation & Filtrage\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#################################################################################################

for i in range(0,len(data)):
    flexion[i][0]=atan2(data[i][2],data[i][1])*180/pi
    flexion[i][1]=atan2(data[i][6],data[i][7])*180/pi
    flexion[i][2]=flexion[i][0]-flexion[i][1]
#if flexion[i][2]>180:
#flexion[i][2]=360-flexion[i][2]
    rotation[i][0]=asin(data[i][0]/sqrt(data[i][0]*data[i][0]+data[i][1]*data[i][1]+data[i][2]*data[i][2]))*180/pi
    rotation[i][1]=asin(-data[i][8]/sqrt(data[i][6]*data[i][6]+data[i][7]*data[i][7]+data[i][8]*data[i][8]))*180/pi
    rotation[i][2]=rotation[i][0]-rotation[i][1]

#plt.plot([flexion[i][2] for i in range(len(flexion))])
flex=[flexion[i][2] for i in range(len(flexion))]
rot=[rotation[i][2] for i in range(len(rotation))]

del rotation

#filtreflex=savgol_filter(flex,67,5)
#filtrerot=savgol_filter(rot,67,5)


filtreflex=savgol_filter(flex,31,2)
filtrerot=savgol_filter(rot,31,2)

##########################################################################################
#//////////////////////////////  Execution    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
##########################################################################################

if exer=="extention":
    filtreflex=savgol_filter(flex,101,5)    
    maximus=extentionAnalyse()
    Series=detectSeries()
    
   
elif exer=="squats_2feet" or exer=="squats_1foot" :
    filtreflex=savgol_filter(flex,101,5)
    exercise="Squat"
    try:       
        maximus=squatAnalyse()
        Series=detectSeries()
        Series=Series[Series[:, 1] != 0]        
    except finish:        
        finish="False"
        print(finish)

elif exer=="proprioception_static" or exer=="proprioception_pillow"  :
    try:
        maximus=Proprio_Analyse()
        Series=detectSeriesProp()    
        Series=Series[Series[:, 1] != 0]        
        if len(Series)==0 and len(maximus)!=0:
            Series=np.zeros([1,2])
            Series[0][0]=maximus[0][0]
            Series[0][1]=maximus[len(maximus)-1][3]
        elif len(Series)==0 and len(maximus)==0:
            Series=np.zeros([1,2])
            Series[0][0]=0
            Series[0][1]=len(flex)-1 
    except finish:
        print("finish")
        
elif exer=="press" :
    filtreflex=savgol_filter(flex,31,2)    
    filtrerot=savgol_filter(rot,31,2)  
    maximus=extentionAnalyse()
    Series=detectSeries()    
    
    
    
    
        
        
        
        
        
   # exportPdf(Series)        
        
   
#################################################################################################
#//////////////////////////////         Mail Service      \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#################################################################################################

Series=np.zeros([2,2])

Series[0][0]=50
Series[0][1]=800
Series[1][0]=900
Series[1][1]=1500
exportPdf(Series)


tmps2=time.clock()
print "Temps d'execution = %f\n" %(tmps2-tmps1)   