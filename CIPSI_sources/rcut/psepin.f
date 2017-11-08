      subroutine psepin
c
c
c     subroutine extraite du programme de toulouse.
c
c
      implicit real*8 (a-h,o-z)
      common/numeriq/xg(64),hg(64)
      dimension fg(13),x(40),h(40),xx(32),hh(32),fgam(13)
      data fg  /  0.1772453850905516d+1,0.1d1,
     1            0.8862269254527581d0,0.1d1,
     2            0.1329340388179137d1,0.2d1,
     3            0.3323350970447843d1,0.6d1,
     4            0.1163172839656745d2,.024d2,
     5            0.5234277778455352d2,0.12d3,
     6            0.2878852778170443d3/
c.......
      data x( 1)/ -.9982377097105589d 00/
      data x( 2)/ -.9907262386994567d 00/
      data x( 3)/ -.9772599499837742d 00/
      data x( 4)/ -.9579168192137913d 00/
      data x( 5)/ -.9328128082786765d 00/
      data x( 6)/ -.9020988069688741d 00/
      data x( 7)/  -.8659595032122592d 00/
      data x( 8)/  -.8246122308333116d 00/
      data x( 9)/  -.7783056514265193d 00/
      data x(10)/  -.7273182551899268d 00/
      data x(11)/  -.6719566846141792d 00/
      data x(12)/  -.6125538896679799d 00/
      data x(13)/  -.5494671250951279d 00/
      data x(14)/  -.4830758016861787d 00/
      data x(15)/  -.4137792043716050d 00/
      data x(16)/  -.3419940908257584d 00/
      data x(17)/  -.2681521850072536d 00/
      data x(18)/  -.1926975807013711d 00/
      data x(19)/  -.1160840706752552d 00/
      data x(20)/  -.3877241750605079d-01/
      data x(21)/   .3877241750605079d-01/
      data x(22)/   .1160840706752552d 00/
      data x(23)/   .1926975807013711d 00/
      data x(24)/   .2681521850072536d 00/
      data x(25)/   .3419940908257584d 00/
      data x(26)/   .4137792043716050d 00/
      data x(27)/   .4830758016861787d 00/
      data x(28)/   .5494671250951279d 00/
      data x(29)/   .6125538896679799d 00/
      data x(30)/   .6719566846141792d 00/
      data x(31)/   .7273182551899268d 00/
      data x(32)/   .7783056514265193d 00/
      data x(33)/   .8246122308333116d 00/
      data x(34)/   .8659595032122592d 00/
      data x(35)/   .9020988069688741d 00/
      data x(36)/   .9328128082786765d 00/
      data x(37)/   .9579168192137913d 00/
      data x(38)/   .9772599499837742d 00/
      data x(39)/   .9907262386994567d 00/
      data x(40)/   .9982377097105589d 00/
c......
      data h( 1)/   .4521277098533191d-02/
      data h( 2)/   .1049828453115281d-01/
      data h( 3)/   .1642105838190788d-01/
      data h( 4)/   .2224584919416695d-01/
      data h( 5)/   .2793700698002340d-01/
      data h( 6)/   .3346019528254784d-01/
      data h( 7)/   .3878216797447202d-01/
      data h( 8)/   .4387090818567327d-01/
      data h( 9)/   .4869580763507223d-01/
      data h(10)/   .5322784698393678d-01/
      data h(11)/   .5743976909939152d-01/
      data h(12)/   .6130624249292891d-01/
      data h(13)/   .6480401345660099d-01/
      data h(14)/   .6791204581523383d-01/
      data h(15)/   .7061164739128672d-01/
      data h(16)/   .7288658239580399d-01/
      data h(17)/   .7472316905796821d-01/
      data h(18)/   .7611036190062619d-01/
      data h(19)/   .7703981816424792d-01/
      data h(20)/   .7750594797842478d-01/
      data h(21)/   .7750594797842478d-01/
      data h(22)/   .7703981816424792d-01/
      data h(23)/   .7611036190062619d-01/
      data h(24)/   .7472316905796821d-01/
      data h(25)/   .7288658239580399d-01/
      data h(26)/   .7061164739128672d-01/
      data h(27)/   .6791204581523383d-01/
      data h(28)/   .6480401345660099d-01/
      data h(29)/   .6130624249292891d-01/
      data h(30)/   .5743976909939152d-01/
      data h(31)/   .5322784698393678d-01/
      data h(32)/   .4869580763507223d-01/
      data h(33)/   .4387090818567327d-01/
      data h(34)/   .3878216797447202d-01/
      data h(35)/   .3346019528254784d-01/
      data h(36)/   .2793700698002340d-01/
      data h(37)/   .2224584919416695d-01/
      data h(38)/   .1642105838190788d-01/
      data h(39)/   .1049828453115281d-01/
      data h(40)/   .4521277098533191d-02/
c......
      data xx( 1)/0.024350292663424430d0/
      data xx( 2)/0.072993121787799040d0/
      data xx( 3)/0.121462819296120600d0/
      data xx( 4)/0.169644420423992800d0/
      data xx( 5)/0.217423643740007100d0/
      data xx( 6)/0.264687162208767400d0/
      data xx( 7)/0.311322871990211000d0/
      data xx( 8)/0.357220158337668100d0/
      data xx( 9)/0.402270157963991600d0/
      data xx(10)/0.446366017253464100d0/
      data xx(11)/0.489403145707053000d0/
      data xx(12)/0.531279464019894500d0/
      data xx(13)/0.571895646202634000d0/
      data xx(14)/0.611155355172393300d0/
      data xx(15)/0.648965471254657300d0/
      data xx(16)/0.685236313054233200d0/
      data xx(17)/0.719881850171610800d0/
      data xx(18)/0.752819907260531900d0/
      data xx(19)/0.783972358943341400d0/
      data xx(20)/0.813265315122797600d0/
      data xx(21)/0.840629296252580400d0/
      data xx(22)/0.865999398154092800d0/
      data xx(23)/0.889315445995114100d0/
      data xx(24)/0.910522137078502800d0/
      data xx(25)/0.929569172131939600d0/
      data xx(26)/0.946411374858402800d0/
      data xx(27)/0.961008799652053700d0/
      data xx(28)/0.973326827789911000d0/
      data xx(29)/0.983336253884626000d0/
      data xx(30)/0.991013371476744300d0/
      data xx(31)/0.996340116771955300d0/
      data xx(32)/0.999305041735772100d0/
c......
      data hh( 1)/0.048690957009139720d0/
      data hh( 2)/0.048575467441503430d0/
      data hh( 3)/0.048344762234802960d0/
      data hh( 4)/0.047999388596458310d0/
      data hh( 5)/0.047540165714830310d0/
      data hh( 6)/0.046968182816210020d0/
      data hh( 7)/0.046284796581314420d0/
      data hh( 8)/0.045491627927418140d0/
      data hh( 9)/0.044590558163756560d0/
      data hh(10)/0.043583724529323450d0/
      data hh(11)/0.042473515123653590d0/
      data hh(12)/0.041262563242623530d0/
      data hh(13)/0.039953741132720340d0/
      data hh(14)/0.038550153178615630d0/
      data hh(15)/0.037055128540240050d0/
      data hh(16)/0.035472213256882380d0/
      data hh(17)/0.033805161837141610d0/
      data hh(18)/0.032057928354851550d0/
      data hh(19)/0.030234657072402480d0/
      data hh(20)/0.028339672614259480d0/
      data hh(21)/0.026377469715054660d0/
      data hh(22)/0.024352702568710870d0/
      data hh(23)/0.022270173808383250d0/
      data hh(24)/0.020134823153530210d0/
      data hh(25)/0.017951715775697340d0/
      data hh(26)/0.015726030476024720d0/
      data hh(27)/0.013463047896718640d0/
      data hh(28)/0.011168139460131130d0/
      data hh(29)/0.008846759826363948d0/
      data hh(30)/0.006504457968978363d0/
      data hh(31)/0.004147033260562468d0/
      data hh(32)/0.001783280721696433d0/
      do 210 i=1,13
  210 fgam(i)=fg(i)
      do 230 i=1,32
      j=i+32
      k=33-i
      xg(i)=-xx(k)
      hg(i)=hh(k)
      xg(j)=xx(i)
  230 hg(j)=hh(i)
      return
      end
