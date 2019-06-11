from contextlib import contextmanager, ExitStack
from itertools import chain
from regex import compile


MAINCHROMS_ENSEMBL = {str(i) for i in range(1, 23)} | {"X", "Y"}
MAINCHROMS_UCSC = {"chr" + s for s in MAINCHROMS_ENSEMBL}
MAINCHROMS_T2T = {"chrX_fixedBionanoSV_centromereV3"}

ALPHABET = list("ACGT")
COMPLEMENTS = dict(zip(ALPHABET, reversed(ALPHABET)))
MOTIF_COMPLEMENTS = {**COMPLEMENTS, **{"[": "]", "]": "[", ".": "."}}
MOTIF_COMPLEMENT_PATTERN = compile(r'|'.join(MOTIF_COMPLEMENTS.keys()))

M19947_HMMGZ = b'+,^C)Jk=m*!X$];:MUq?\'ZBJ;.Jb%RHM[U-p*$Z&-/48WE[8U1KH(,<N"1"RE=R_Npe7m(M#o"t6-.OmULs77[T;[B7K$?%kPP5BPs52Xr2<iQf()?,qY0`&MG1GV[(4iV:J9BmhCU\\LV\\dZOOc[h:pVV8*<pO@Arpk=`:G6m\\X0XQgmsj$WpR`3KQP\\8MI!kP?k2cC6h!SaaJ*?.bB"mliYBOQ-qm"X%<M,(3d6:mID_A1g^Fc:":=%3<k-*P(E,A@a\\c-JRe@N(BPKidFmcM*N?_6W0bPq?<3F%WFD)HX*12V))[ssXUH?""Pp!nstk1c,:T.`:<]_[)$olSu,ghl>+QT`*5L\\KEg\\2ZK$SfJV+oD_O1E4pVjjh\'R4m^RLi=POph\\iW+KqY8.EGWjsHjDt:kZAbE0b;AHcFK9SHhn,8:7s+`)o`qZ\\5LMmp^FUE?:YR>BXVBSpGWoB^GnFedB)_aE*/_dZr8\'(jh5;<#[+D;/Z,e8<0@cO4^--55G-(X=*fY[LMr8jkobM[4@[#-mq>ME:#SmW\\$4\'Z:_03>L8"ItLOmUZ"*F<M6D21j?HB.`??[KE!@))^JVaM$2i[q&qc.m@3n.n.)O0roX(E@t`?kAiYi+*DSJ4%p*%UB;BG65U\'>2CRfg@`\'GXipJ.>eBtVXo\\kE#/LWP<X-CKr6W\\q#S*C-Jt52H!*pcc\\cgr5Q,2^"!sCOe_&X*B#mAGYr6W]rlQ8.g17^+03g&N0@]j"O9*jHt_VfNFJ7)ju$tMD-aE/\'1TaQTBQ[n")&DVdgd29Fti,lioF?1ui4P(2O"u;`)A3\'nXV)$Cj"mDd<cq=Q+V)$Cj"e8=*S5\'*2W=8DsJP;sMR8*d/W=8Ds!IC^Ro[3_cBajX^!2hUWV_`j0)?[Z0+O_2qVDEa/)?[Z0J6j,r]ek7G)?[Z0_!-fokmX5O`]u0daE5F`ko?@_c9O$\'TU3qh0--^Pc9O$\'E&qWN3X=lF1XJc"17^433X=lF1XJc"kbU)^9*a[W1XJc"SB`"B/c^28%.l<*9G\\,*/c^28%.l:TH5R.`1&uV<%.l:TV[E^,$\\2Of"2,*cdfi+9$\\2Of"2,+n?3UMF9RpFr"2.CDl&rmI&tkT55VEelFNKX/&tkT55VEeDQcVKG\';1]65VEdi\\n#ck"R\'DPd)uPI>r#s%"uomuE4Q65?"7ib#%Y<bOQi28>=\'49"@+UZ"2/M9VJ?U.XoXrn0k@j"d*\'Eo0,p#;)FLA@/R`XckQG1c"8&lTd)unY?"\\,f%V3/rOQeek0099Llp.H6%%49qnh)ZK&p+GqkrOIGFWf.UL.p3M9\\1.;>mKS)-_QfK>5jAC;+i%A0(\'@Ip3(l#bFNP#^0P*M!eQ\\\'YQ9?86*LHHc5@*t!.ohhDus;h)7i3<p2FNo_]9:^5VHnYT\'8]rJ`,X^]c=HG+1`qP?ssVZS6h60ot*TqS?=V,?$/E:1C+IUE%hPmi6WdfTkf2:`ZVbn?%RT5<WS^c(n=>N+Lj04"u;`\'A3*jg01YF+LJA6hC^bo^e&tT]cm8[i\'*b@u"9^*)?m*o)6`-%W*\'+:MW;sV2#(!6K]i52^WIWGHK@;8709q)p.ZinZ+k%l*09q)p.R5YPJR0fk09q)p.c>Eki,f=Q94ac_#%Tc)aE.X594OW]#1-"QTU,l094OW]"ur/.31U*BdMV2E!L_4M17\\I<dMV2E!VFeWkbS=PdMV2E"&A"f:?MVjF9[\'U5c/UW9]lDhF9[\'U!53IYHKbFCF9[\'UJC?N2dlBZHSUk*.OQi8YdlBZHSUk)O5VIb@?3WYJK7R>5i,lfnF90b3@%`heA@c0iF90b3@%`hed*#giQN>IV@%`he3*e3t>laLP)!GN4R4]@4>lsXR)<bW5oJ.<JAHMKZ)<bW5:6uAaQN4]m#\'u\\&VFuY;QN4]m#\'u[;]hAWVR/joo#("r&VeZ]9#L#`<$o:dr"#=.laE.laQ]6Q&2Lo*l#("r&HJ\'/!,qpR?1<`:g%A9&IJR15WdkukRf)U-M)!GN9k[eD"QQ;`cSUk)K1K#^?N0N#sHJ\'2"AHstrB?508TU24#??Uj*Fp<9G!IEE:kTqo-V)$7cW?%XuSUb%hZolfs=ouO`l/elBYeFKG":Kk?jn*TChQANaJ3@]WHe8A8c6,(H_8.6RJO>b4LTUZ0)h[&YF9r"p.DUlV@@30$[/tq1K4"]_$\\RC!oQ&HaJ*<m5JrN(cEdB`\\$+^8tF:9s\'Ua@3D)$Dl`E+9I7JO;pm01D2Ub>CnTGWEGW[Hi#C?q-#-bX`=C5d3!L"u;`\'A3+_gcoK";$UAYb#Yq$c,X[muS;q#IMF[sT1XJbu17^123X+`D1=/YtkbU&]9*OOU1=/[%7Z&nD8.c$hV)$Cj"YuW(A3\'nXV)$Cj"mDd<cq=Q+V)$Cj?q-#-bX`<JkbU#\\9*6r?0k<<IkbU#\\9415H0[*\'""#?E%A3\'fJSM=Y9"#?E%A3\'fJSM=Y9"#?E%A3+*fMD37F&I1(^1LrNO17^433X7F21LrNOSB`"Bc@J[P%%.TmOWk2i;T+,[k7!d(JP>58REg@Wk7!b2!2g2/V[FWr"[f-hAu*Lt1*]5!r)Mp];t]WWi5>6,k60,Mj8,fgoD%k<:N129"`_!O!!'


@contextmanager
def ReadFileChain(filenames, manager):
    """Chain records from all filenames in list bams, replicating behavior of pysam context managers"""
    with ExitStack() as stack:
        yield chain(*(
            stack.enter_context(manager(filename))
            for filename in filenames
        ))


def motif_revcomp(motif, ignorecase=True):
    """Reverse-complement a regex-like motif; only allows a subset of regex syntax (ACGT, dot, square brackets)"""
    if ignorecase:
        matcher = lambda match: MOTIF_COMPLEMENTS[match.group().upper()]
    else:
        matcher = lambda match: MOTIF_COMPLEMENTS[match.group()]
    try:
        return MOTIF_COMPLEMENT_PATTERN.sub(matcher, motif[::-1])
    except KeyError:
        raise ValueError("Unsupported character(s) in motif: {}".format(motif))
