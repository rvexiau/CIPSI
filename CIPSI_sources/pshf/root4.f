      subroutine root4
      implicit real*8 (a-h,o-z)
c          *****   version february 16,1975   *****
      common/root/x,u(9),w(9),nroots
      equivalence (u(1),rt1),(u(2),rt2),(u(3),rt3),(u(4),rt4),(u(5),rt5)
      equivalence (w(1),ww1),(w(2),ww2),(w(3),ww3),(w(4),ww4),(w(5),ww5)
      data r14,pie4/1.45303521503316d-01, 7.85398163397448d-01/
      data r24,w24/ 1.33909728812636d+00, 2.34479815323517d-01/
      data r34,w34/ 3.92696350135829d+00, 1.92704402415764d-02/
      data r44,w44/ 8.58863568901199d+00, 2.25229076750736d-04/
      if(x.gt.15.0d+00) go to 470
      if(x.gt.5.0d+00) go to 450
      if(x.gt.1.0d+00) go to 430
      if(x.gt.3.0d-7) go to 420
c     x is approximately zero.                   nroots = 4
      rt1=3.48198973061471d-02    -4.09645850660395d-03 *x
      rt2=3.81567185080042d-01    -4.48902570656719d-02 *x
      rt3=1.73730726945891d+00    -2.04389090547327d-01 *x
      rt4=1.18463056481549d+01    -1.39368301742312d+00 *x
      ww1=3.62683783378362d-01    -3.13844305713928d-02 *x
      ww2=3.13706645877886d-01    -8.98046242557724d-02 *x
      ww3=2.22381034453372d-01    -1.29314370958973d-01 *x
      ww4=1.01228536290376d-01    -8.28299075414321d-02 *x
      return
c     x=0.0 to 1.0                               nroots = 4
  420 rt1=            ((((((-1.95309614628539d-10*x+5.19765728707592d-09
     1)*x-1.01756452250573d-07 )*x+1.72365935872131d-06
     2)*x-2.61203523522184d-05 )*x+3.52921308769880d-04
     3)*x-4.09645850658433d-03 )*x+3.48198973061469d-02
      rt2=             (((((-1.89554881382342d-08*x+3.07583114342365d-07
     1)*x+1.270981734393d-06)*x-1.417298563884d-04)*x+3.226979163176d-03
     2)*x-4.48902570678178d-02 )*x+3.81567185080039d-01
      rt3=            (((((( 1.77280535300416d-09*x+3.36524958870615d-08
     1)*x-2.58341529013893d-07 )*x-1.13644895662320d-05
     2)*x-7.91549618884063d-05 )*x+1.03825827346828d-02
     3)*x-2.04389090525137d-01 )*x+1.73730726945889d+00
      rt4=             (((((-5.61188882415248d-08*x-2.49480733072460d-07
     1)*x+3.428685057114d-06)*x+1.679007454539d-04)*x+4.722855585715d-02
     2)*x-1.39368301737828d+00 )*x+1.18463056481543d+01
      ww1=            ((((((-1.14649303201279d-08*x+1.88015570196787d-07
     1)*x-2.33305875372323d-06 )*x+2.68880044371597d-05
     2)*x-2.94268428977387d-04 )*x+3.06548909776613d-03
     3)*x-3.13844305680096d-02 )*x+3.62683783378335d-01
      ww2=          ((((((((-4.11720483772634d-09*x+6.54963481852134d-08
     1)*x-7.20045285129626d-07 )*x+6.93779646721723d-06
     2)*x-6.05367572016373d-05 )*x+4.74241566251899d-04
     3)*x-3.26956188125316d-03 )*x+1.91883866626681d-02
     4)*x-8.98046242565811d-02 )*x+3.13706645877886d-01
      ww3=          ((((((((-3.41688436990215d-08*x+5.07238960340773d-07
     1)*x-5.01675628408220d-06 )*x+4.20363420922845d-05
     2)*x-3.08040221166823d-04 )*x+1.94431864731239d-03
     3)*x-1.02477820460278d-02 )*x+4.28670143840073d-02
     4)*x-1.29314370962569d-01 )*x+2.22381034453369d-01
      ww4=         ((((((((( 4.99660550769508d-09*x-7.94585963310120d-08
     1)*x+8.359072409485d-07)*x-7.422369210610d-06)*x+5.763374308160d-05
     2)*x-3.86645606718233d-04 )*x+2.18417516259781d-03
     3)*x-9.99791027771119d-03 )*x+3.48791097377370d-02
     4)*x-8.28299075413889d-02 )*x+1.01228536290376d-01
      return
c     x= 1.0 to 5.0                              nroots = 4
  430 y=x-3.0d+00
      rt1=         (((((((((-1.48570633747284d-15*y-1.33273068108777d-13
     1)*y+4.068543696670d-12)*y-9.163164161821d-11)*y+2.046819017845d-09
     2)*y-4.03076426299031d-08 )*y+7.29407420660149d-07
     3)*y-1.23118059980833d-05 )*y+1.88796581246938d-04
     4)*y-2.53262912046853d-03 )*y+2.51198234505021d-02
      rt2=         ((((((((( 1.35830583483312d-13*y-2.29772605964836d-12
     1)*y-3.821500128045d-12)*y+6.844424214735d-10)*y-1.048063352259d-08
     2)*y+1.50083186233363d-08 )*y+3.48848942324454d-06
     3)*y-1.08694174399193d-04 )*y+2.08048885251999d-03
     4)*y-2.91205805373793d-02 )*y+2.72276489515713d-01
      rt3=         ((((((((( 5.02799392850289d-13*y+1.07461812944084d-11
     1)*y-1.482277886411d-10)*y-2.153585661215d-09)*y+3.654087802817d-08
     2)*y+5.15929575830120d-07 )*y-9.52388379435709d-06
     3)*y-2.16552440036426d-04 )*y+9.03551469568320d-03
     4)*y-1.45505469175613d-01 )*y+1.21449092319186d+00
      rt4=         (((((((((-1.08510370291979d-12*y+6.41492397277798d-11
     1)*y+7.542387436125d-10)*y-2.213111836647d-09)*y-1.448228963549d-07
     2)*y-1.95670833237101d-06 )*y-1.07481314670844d-05
     3)*y+1.49335941252765d-04 )*y+4.87791531990593d-02
     4)*y-1.10559909038653d+00 )*y+8.09502028611780d+00
      ww1=        ((((((((((-4.65801912689961d-14*y+7.58669507106800d-13
     1)*y-1.186387548048d-11)*y+1.862334710665d-10)*y-2.799399389539d-09
     2)*y+4.148972684255d-08)*y-5.933568079600d-07)*y+8.168349266115d-06
     3)*y-1.08989176177409d-04 )*y+1.41357961729531d-03
     4)*y-1.87588361833659d-02 )*y+2.89898651436026d-01
      ww2=      ((((((((((((-1.46345073267549d-14*y+2.25644205432182d-13
     1)*y-3.116258693847d-12)*y+4.321908756610d-11)*y-5.673270062669d-10
     2)*y+7.006295962960d-09)*y-8.120186517000d-08)*y+8.775294645770d-07
     3)*y-8.77829235749024d-06 )*y+8.04372147732379d-05
     4)*y-6.64149238804153d-04 )*y+4.81181506827225d-03
     5)*y-2.88982669486183d-02 )*y+1.56247249979288d-01
      ww3=     ((((((((((((( 9.06812118895365d-15*y-1.40541322766087d-13
     1)*y+1.919270015269d-12)*y-2.605135739010d-11)*y+3.299685839012d-10
     2)*y-3.86354139348735d-09 )*y+4.16265847927498d-08
     3)*y-4.09462835471470d-07 )*y+3.64018881086111d-06
     4)*y-2.88665153269386d-05 )*y+2.00515819789028d-04
     5)*y-1.18791896897934d-03 )*y+5.75223633388589d-03
     6)*y-2.09400418772687d-02 )*y+4.85368861938873d-02
      ww4=    ((((((((((((((-9.74835552342257d-16*y+1.57857099317175d-14
     1)*y-2.249993780112d-13)*y+3.173422008953d-12)*y-4.161159459680d-11
     2)*y+5.021343560166d-10)*y-5.545047534808d-09)*y+5.554146993491d-08
     3)*y-4.99048696190133d-07 )*y+3.96650392371311d-06
     4)*y-2.73816413291214d-05 )*y+1.60106988333186d-04
     5)*y-7.64560567879592d-04 )*y+2.81330044426892d-03
     6)*y-7.16227030134947d-03 )*y+9.66077262223353d-03
      return
  450 if(x.gt.10.0d+00) go to 460
c     x=5.0 to 10.0                              nroots = 4
      y=x-7.5d+00
      rt1=         ((((((((( 4.64217329776215d-15*y-6.27892383644164d-15
     1)*y+3.462236347446d-13)*y-2.927229355350d-11)*y+5.090355371676d-10
     2)*y-9.97272656345253d-09 )*y+2.37835295639281d-07
     3)*y-4.60301761310921d-06 )*y+8.42824204233222d-05
     4)*y-1.37983082233081d-03 )*y+1.66630865869375d-02
      rt2=         ((((((((( 2.93981127919047d-14*y+8.47635639065744d-13
     1)*y-1.446314544774d-11)*y-6.149155555753d-12)*y+8.484275604612d-10
     2)*y-6.10898827887652d-08 )*y+2.39156093611106d-06
     3)*y-5.35837089462592d-05 )*y+1.00967602595557d-03
     4)*y-1.57769317127372d-02 )*y+1.74853819464285d-01
      rt3=        (((((((((( 2.93523563363000d-14*y-6.40041776667020d-14
     1)*y-2.695740446312d-12)*y+1.027082960169d-10)*y-5.822038656780d-10
     2)*y-3.159991002539d-08)*y+4.327249251331d-07)*y+4.856768455119d-06
     3)*y-2.54617989427762d-04 )*y+5.54843378106589d-03
     4)*y-7.95013029486684d-02 )*y+7.20206142703162d-01
      rt4=       (((((((((((-1.62212382394553d-14*y+7.68943641360593d-13
     1)*y+5.764015756615d-12)*y-1.380635298784d-10)*y-1.476849808675d-09
     2)*y+1.84347052385605d-08 )*y+3.34382940759405d-07
     3)*y-1.39428366421645d-06 )*y-7.50249313713996d-05
     4)*y-6.26495899187507d-04 )*y+4.69716410901162d-02
     5)*y-6.66871297428209d-01 )*y+4.11207530217806d+00
      ww1=        ((((((((((-1.65995045235997d-15*y+6.91838935879598d-14
     1)*y-9.131223418888d-13)*y+1.403341829454d-11)*y-3.672235069444d-10
     2)*y+6.366962546990d-09)*y-1.039220021671d-07)*y+1.959098751715d-06
     3)*y-3.33474893152939d-05 )*y+5.72164211151013d-04
     4)*y-1.05583210553392d-02 )*y+2.26696066029591d-01
      ww2=      ((((((((((((-3.57248951192047d-16*y+6.25708409149331d-15
     1)*y-9.657033089714d-14)*y+1.507864898748d-12)*y-2.332522256110d-11
     2)*y+3.428545616603d-10)*y-4.698730937661d-09)*y+6.219977635130d-08
     3)*y-7.83008889613661d-07 )*y+9.08621687041567d-06
     4)*y-9.86368311253873d-05 )*y+9.69632496710088d-04
     5)*y-8.14594214284187d-03 )*y+8.50218447733457d-02
      ww3=     ((((((((((((( 1.64742458534277d-16*y-2.68512265928410d-15
     1)*y+3.788890667676d-14)*y-5.508918529823d-13)*y+7.555896810069d-12
     2)*y-9.69039768312637d-11 )*y+1.16034263529672d-09
     3)*y-1.28771698573873d-08 )*y+1.31949431805798d-07
     4)*y-1.23673915616005d-06 )*y+1.04189803544936d-05
     5)*y-7.79566003744742d-05 )*y+5.03162624754434d-04
     6)*y-2.55138844587555d-03 )*y+1.13250730954014d-02
      ww4=    ((((((((((((((-1.55714130075679d-17*y+2.57193722698891d-16
     1)*y-3.626606654097d-15)*y+5.234734676175d-14)*y-7.067105402134d-13
     2)*y+8.793512664890d-12)*y-1.006088923498d-10)*y+1.050565098393d-09
     3)*y-9.91517881772662d-09 )*y+8.35835975882941d-08
     4)*y-6.19785782240693d-07 )*y+3.95841149373135d-06
     5)*y-2.11366761402403d-05 )*y+9.00474771229507d-05
     6)*y-2.78777909813289d-04 )*y+5.26543779837487d-04
      return
c     x=10.0 to 15.0                             nroots = 4
  460 y=x-12.5d+00
      rt1=       ((((((((((( 4.94869622744119d-17*y+8.03568805739160d-16
     1)*y-5.599125915431d-15)*y-1.378685560217d-13)*y+7.006511663249d-13
     2)*y+1.30391406991118d-11 )*y+8.06987313467541d-11
     3)*y-5.20644072732933d-09 )*y+7.72794187755457d-08
     4)*y-1.61512612564194d-06 )*y+4.15083811185831d-05
     5)*y-7.87855975560199d-04 )*y+1.14189319050009d-02
      rt2=       ((((((((((( 4.89224285522336d-16*y+1.06390248099712d-14
     1)*y-5.446260182933d-14)*y-1.613630106295d-12)*y+3.910179118937d-12
     2)*y+1.90712434258806d-10 )*y+8.78470199094761d-10
     3)*y-5.97332993206797d-08 )*y+9.25750831481589d-07
     4)*y-2.02362185197088d-05 )*y+4.92341968336776d-04
     5)*y-8.68438439874703d-03 )*y+1.15825965127958d-01
      rt3=        (((((((((( 6.12419396208408d-14*y+1.12328861406073d-13
     1)*y-9.051094103059d-12)*y-4.781797525341d-11)*y+1.660828868694d-09
     2)*y+4.499058798868d-10)*y-2.519549641933d-07)*y+4.977444040180d-06
     3)*y-1.25858350034589d-04 )*y+2.70279176970044d-03
     4)*y-3.99327850801083d-02 )*y+4.33467200855434d-01
      rt4=       ((((((((((( 4.63414725924048d-14*y-4.72757262693062d-14
     1)*y-1.001926833832d-11)*y+6.074107718414d-11)*y+1.576976911942d-09
     2)*y-2.01186401974027d-08 )*y-1.84530195217118d-07
     3)*y+5.02333087806827d-06 )*y+9.66961790843006d-06
     4)*y-1.58522208889528d-03 )*y+2.80539673938339d-02
     5)*y-2.78953904330072d-01 )*y+1.82835655238235d+00
      ww4=     ((((((((((((( 2.90401781000996d-18*y-4.63389683098251d-17
     1)*y+6.274018198326d-16)*y-8.936002188168d-15)*y+1.194719074934d-13
     2)*y-1.45501321259466d-12 )*y+1.64090830181013d-11
     3)*y-1.71987745310181d-10 )*y+1.63738403295718d-09
     4)*y-1.39237504892842d-08 )*y+1.06527318142151d-07
     5)*y-7.27634957230524d-07 )*y+4.12159381310339d-06
     6)*y-1.74648169719173d-05 )*y+8.50290130067818d-05
      ww3=      ((((((((((((-4.19569145459480d-17*y+5.94344180261644d-16
     1)*y-1.148797566469d-14)*y+1.881303962576d-13)*y-2.413554618391d-12
     2)*y+3.372127423047d-11)*y-4.933988617784d-10)*y+6.116545396281d-09
     3)*y-6.69965691739299d-08 )*y+7.52380085447161d-07
     4)*y-8.08708393262321d-06 )*y+6.88603417296672d-05
     5)*y-4.67067112993427d-04 )*y+5.42313365864597d-03
      ww2=        ((((((((((-6.22272689880615d-15*y+1.04126809657554d-13
     1)*y-6.842418230913d-13)*y+1.576841731919d-11)*y-4.203948834175d-10
     2)*y+6.287255934781d-09)*y-8.307159819228d-08)*y+1.356478091922d-06
     3)*y-2.08065576105639d-05 )*y+2.52396730332340d-04
     4)*y-2.94484050194539d-03 )*y+6.01396183129168d-02
      ww1=       (((-1.8784686463512d-01/x+2.2991849164985d-01)/x
     1-4.9893752514047d-01)/x-2.1916512131607d-05)*exp(-x)
     2 +sqrt(pie4/x)-ww4-ww3-ww2
      return
  470 ww1=sqrt(pie4/x)
      if(x.gt.35.0d+00) go to 490
      if(x.gt.20.0d+00) go to 480
c     x=15.0 to 20.0                             nroots = 4
      y=x-17.5d+00
      rt1=       ((((((((((( 4.36701759531398d-17*y-1.12860600219889d-16
     1)*y-6.149849164164d-15)*y+5.820231579541d-14)*y+4.396602872143d-13
     2)*y-1.24330365320172d-11 )*y+6.71083474044549d-11
     3)*y+2.43865205376067d-10 )*y+1.67559587099969d-08
     4)*y-9.32738632357572d-07 )*y+2.39030487004977d-05
     5)*y-4.68648206591515d-04 )*y+8.34977776583956d-03
      rt2=       ((((((((((( 4.98913142288158d-16*y-2.60732537093612d-16
     1)*y-7.775156445127d-14)*y+5.766105220086d-13)*y+6.432696729600d-12
     2)*y-1.39571683725792d-10 )*y+5.95451479522191d-10
     3)*y+2.42471442836205d-09 )*y+2.47485710143120d-07
     4)*y-1.14710398652091d-05 )*y+2.71252453754519d-04
     5)*y-4.96812745851408d-03 )*y+8.26020602026780d-02
      rt3=       ((((((((((( 1.91498302509009d-15*y+1.48840394311115d-14
     1)*y-4.316925145767d-13)*y+1.186495793471d-12)*y+4.615806713055d-11
     2)*y-5.54336148667141d-10 )*y+3.48789978951367d-10
     3)*y-2.79188977451042d-09 )*y+2.09563208958551d-06
     4)*y-6.76512715080324d-05 )*y+1.32129867629062d-03
     5)*y-2.05062147771513d-02 )*y+2.88068671894324d-01
      rt4=       (((((((((((-5.43697691672942d-15*y-1.12483395714468d-13
     1)*y+2.826607936174d-12)*y-1.266734493280d-11)*y-4.258722866437d-10
     2)*y+9.45486578503261d-09 )*y-5.86635622821309d-08
     3)*y-1.28835028104639d-06 )*y+4.41413815691885d-05
     4)*y-7.61738385590776d-04 )*y+9.66090902985550d-03
     5)*y-1.01410568057649d-01 )*y+9.54714798156712d-01
      ww4=      ((((((((((((-7.56882223582704d-19*y+7.53541779268175d-18
     1)*y-1.157318032236d-16)*y+2.411195002314d-15)*y-3.601794386996d-14
     2)*y+4.082150659615d-13)*y-4.289542980767d-12)*y+5.086829642731d-11
     3)*y-6.35435561050807d-10 )*y+6.82309323251123d-09
     4)*y-5.63374555753167d-08 )*y+3.57005361100431d-07
     5)*y-2.40050045173721d-06 )*y+4.94171300536397d-05
      ww3=       (((((((((((-5.54451040921657d-17*y+2.68748367250999d-16
     1)*y+1.349020069254d-14)*y-2.507452792892d-13)*y+1.944339743818d-12
     2)*y-1.29816917658823d-11 )*y+3.49977768819641d-10
     3)*y-8.67270669346398d-09 )*y+1.31381116840118d-07
     4)*y-1.36790720600822d-06 )*y+1.19210697673160d-05
     5)*y-1.42181943986587d-04 )*y+4.12615396191829d-03
      ww2=       (((((((((((-1.86506057729700d-16*y+1.16661114435809d-15
     1)*y+2.563712856363d-14)*y-4.498350984631d-13)*y+1.765194089338d-12
     2)*y+9.04483676345625d-12 )*y+4.98930345609785d-10
     3)*y-2.11964170928181d-08 )*y+3.98295476005614d-07
     4)*y-5.49390160829409d-06 )*y+7.74065155353262d-05
     5)*y-1.48201933009105d-03 )*y+4.97836392625268d-02
      ww1=        (( 1.9623264149430d-01/x-4.9695241464490d-01)/x
     1-6.0156581186481d-05)*exp(-x)+ww1-ww2-ww3-ww4
      return
c     x=20.0 to 35.0                             nroots = 4
  480 e=exp(-x)
      rt1=    ((((((-4.45711399441838d-05*x+1.27267770241379d-03)*x
     1 -2.36954961381262d-01)*x+1.54330657903756d+01)*x
     2 -5.22799159267808d+02)*x+1.05951216669313d+04)*x
     3 +       (-2.51177235556236d+06/x+8.72975373557709d+05)/x
     4 -1.29194382386499d+05)*e + r14/(x-r14)
      rt2=     (((((-7.85617372254488d-02*x+6.35653573484868d+00)*x
     1 -3.38296938763990d+02)*x+1.25120495802096d+04)*x
     2 -3.16847570511637d+05)*x +        ((-1.02427466127427d+09/x
     3 +3.70104713293016d+08)/x-5.87119005093822d+07)/x
     4 +5.38614211391604d+06)*e + r24/(x-r24)
      rt3=     (((((-2.37900485051067d-01*x+1.84122184400896d+01)*x
     1 -1.00200731304146d+03)*x+3.75151841595736d+04)*x
     2 -9.50626663390130d+05)*x +        ((-2.88139014651985d+09/x
     3 +1.06625915044526d+09)/x-1.72465289687396d+08)/x
     4 +1.60419390230055d+07)*e + r34/(x-r34)
      rt4=    ((((((-6.00691586407385d-04*x-3.64479545338439d-01)*x
     1 +1.57496131755179d+01)*x-6.54944248734901d+02)*x
     2 +1.70830039597097d+04)*x-2.90517939780207d+05)*x
     3 +       (+3.49059698304732d+07/x-1.64944522586065d+07)/x
     4 +2.96817940164703d+06)*e + r44/(x-r44)
      if(x.le.25.0d+00)
     1ww4=   ((((((( 2.33766206773151d-07*x-3.81542906607063d-05)*x
     1 +3.51416601267000d-03)*x-1.66538571864728d-01)*x
     2 +4.80006136831847d+00)*x-8.73165934223603d+01)*x
     3 +9.77683627474638d+02)*x +           1.66000945117640d+04/x
     4 -6.14479071209961d+03)*e + w44*ww1
      if(x.gt.25.0d+00)
     1ww4=    (((((( 5.74245945342286d-06*x-7.58735928102351d-05)*x
     1 +2.35072857922892d-04)*x-3.78812134013125d-03)*x
     2 +3.09871652785805d-01)*x-7.11108633061306d+00)*x
     3 +5.55297573149528d+01)*e + w44*ww1
      ww3=    (((((( 2.36392855180768d-04*x-9.16785337967013d-03)*x
     1 +4.62186525041313d-01)*x-1.96943786006540d+01)*x
     2 +4.99169195295559d+02)*x-6.21419845845090d+03)*x
     3 +      ((+5.21445053212414d+07/x-1.34113464389309d+07)/x
     4 +1.13673298305631d+06)/x-2.81501182042707d+03)*e + w34*ww1
      ww2=    (((((( 7.29841848989391d-04*x-3.53899555749875d-02)*x
     1 +2.07797425718513d+00)*x-1.00464709786287d+02)*x
     2 +3.15206108877819d+03)*x-6.27054715090012d+04)*x
     3 +       (+1.54721246264919d+07/x-5.26074391316381d+06)/x
     4 +7.67135400969617d+05)*e + w24*ww1
      ww1=        (( 1.9623264149430d-01/x-4.9695241464490d-01)/x
     1-6.0156581186481d-05)*e + ww1-ww2-ww3-ww4
      return
  490 if(x.gt.53.0d+00) go to 495
c     x=35.0 to 53.0                             nroots = 4
      e=exp(-x)*(x*x)**2
      rt4=        ((-2.19135070169653d-03*x-1.19108256987623d-01)*x
     1 -7.50238795695573d-01)*e + r44/(x-r44)
      rt3=        ((-9.65842534508637d-04*x-4.49822013469279d-02)*x
     1 +6.08784033347757d-01)*e + r34/(x-r34)
      rt2=        ((-3.62569791162153d-04*x-9.09231717268466d-03)*x
     1 +1.84336760556262d-01)*e + r24/(x-r24)
      rt1=        ((-4.07557525914600d-05*x-6.88846864931685d-04)*x
     1 +1.74725309199384d-02)*e + r14/(x-r14)
      ww4=        (( 5.76631982000990d-06*x-7.89187283804890d-05)*x
     1 +3.28297971853126d-04)*e + w44*ww1
      ww3=        (( 2.08294969857230d-04*x-3.77489954837361d-03)*x
     1 +2.09857151617436d-02)*e + w34*ww1
      ww2=        (( 6.16374517326469d-04*x-1.26711744680092d-02)*x
     1 +8.14504890732155d-02)*e + w24*ww1
      ww1=ww1-ww2-ww3-ww4
      return
c     x=47.0 to infinity                         nroots = 4
  495 rt1=r14/(x-r14)
      rt2=r24/(x-r24)
      rt3=r34/(x-r34)
      rt4=r44/(x-r44)
      ww4=w44*ww1
      ww3=w34*ww1
      ww2=w24*ww1
      ww1=ww1-ww2-ww3-ww4
      return
      end
