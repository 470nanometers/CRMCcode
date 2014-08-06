#!/usr/bin/python
#converting .py to .pyc
### This should read in a PDG particle ID, and write out a CORSIKA ID for the particle

def convert(ID):  #gonna be a lot of if-> then,  also need specifics for heavy ions, I need CORSIKA IDs through number 69 I think, beyond that are decays and charm/bottom
  f = open('/home/dlahurd/Desktop/nomatch-B.txt', 'a+')
  inid = int(ID)
  #photons
  if inid == 22:
    corID = str(1)
  #positron
  elif inid == -11:
    corID = str(2)
  #electron
  elif inid == 11:
    corID = str(3)
  #NONE
  #elif inid == :
    #corID = str(4)
  #mu+
  elif inid == -13:
    corID = str(5)
  #mu-
  elif inid == 13:
    corID = str(6)
  #pi0
  elif inid == 111:
    corID = str(7)
  #pi+
  elif inid == 211:
    corID = str(8)
  #pi-
  elif inid == -211:
    corID = str(9)
  #K0L
  elif inid == 130:
    #print 'K0L'
    corID = str(10)
  #K+
  elif inid == 321:
    #print 'K+'
    corID = str(11)
  #K-
  elif inid == -321:
    #print 'K-'
    corID = str(12)
  #neutron
  elif inid == 2112:
    corID = str(13)
  #proton
  elif inid == 2212:
    corID = str(14)
  #anti-proton
  elif inid == -2212:
    corID = str(15)
  #K0S
  elif inid == 310:
    #print 'K0S'
    corID = str(16)
  #eta
  elif inid == 221:
    corID = str(17)
  #Lambda
  elif inid == 3122:
    corID = str(18)
  #Sigma+
  elif inid == 3222:
    corID = str(19)
  #Sigma0
  elif inid == 3212:
    corID = str(20)
  #Sigma-
  elif inid == 3112:
    corID = str(21)
  #Xi0
  elif inid == 3322:
    corID = str(22)
  #Xi-
  elif inid == 3312:
    corID = str(23)
  #Omega-
  elif inid == 3334:
    corID = str(24)
  #anti-neutron
  elif inid == -2112:
    corID = str(25)
  #anti-Lambda
  elif inid == -3122:
    corID = str(26)
  #anti-Sigma-
  elif inid == -3222: #assuming this is the counter part to Sigma+
    corID = str(27)
  #anti-Sigma0
  elif inid == -3212:
    corID = str(28)
  #anti-Sigma+
  elif inid == -3112:
    corID = str(29)
  #anti-Xi0
  elif inid == -3322:
    corID = str(30)
  #anti-Xi+
  elif inid == -3312:
    corID = str(31)
  #anti-Omega+
  elif inid == -3334:
    corID = str(32)
  #omega
  elif inid == 223: #omega782
    corID = str(50)
  #rho0
  elif inid == 113: #rho770
    corID = str(51)
  #rho+
  elif inid == 213:
    corID = str(52)
  #rho-
  elif inid == -213:
    corID = str(53)
  #Delta++
  elif inid == 2224:
    corID = str(54)
  #Delta+
  elif inid == 2214:
    corID = str(55)
  #Delta0
  elif inid == 2114:
    corID = str(56)
  #Delta-
  elif inid == 1114:
    corID = str(57)
  #anti-Delta--
  elif inid == -2224:
    corID = str(58)
  #anti-Delta-
  elif inid == -2214:
    corID = str(59)
  #anti-Delta0
  elif inid == -2114:
    corID = str(60)
  #anti-Delta+
  elif inid == -1114:
    corID = str(61)
  #K*0
  elif inid == 313:
    corID = str(62)
  #K*+
  elif inid == 323:
    corID = str(63)
  #K*-
  elif inid == -323:
    corID = str(64)
  #anti-K*0
  elif inid == -313:
    corID = str(65)
  #electron neutrino
  elif inid == 12:
    corID = str(66)
  # anti-e neutron0
  elif inid == -12:
    corID = str(67)
  #muon neutrino
  elif inid == 14:
    corID = str(68)
  #anti-mu neutrino
  elif inid == -14:
    corID = str(69)

  elif inid > 1000000000:
    #print inid
    z = (inid -1000000000)/10000 #need the 10000ths digits, since this is integer division i think this will work
    a = (inid -1000000000 -(z*10000))/10 
    corID = (a * 100) + z

  else: 
    corID =-1 #if no match return -1
    #print 'NO MATCH', inid
    f.write(str(inid)+'\n')

  f.close()
  return corID

### NO CORSIKA EQUIVALENT FOR:
# 12122 is Delta 1900
# 3224 is Sigma*+ aka Sigma 1385, 
# 3214 is Sigma*0
# 333 is phi
# 311 is K0
# 3314 is Xi*-
# 1216 is ???
# 1214 is ???
# 22124 is ???
# 33122 is
# 331 is
# 21212 is
# 3114 is
# 1112 is
# 1116 is
# 12224 is
# 3314 is 
# 12126 is 
# 21214 is 
# 22214 3114 3324 31214 31114 11114 1212 23212 22212 12112 12126 22114 13116 22122 23112 32124 12212
# 
