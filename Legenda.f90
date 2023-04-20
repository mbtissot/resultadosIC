!Decay zerado:
	Funcx_RcdLLS= ZL(Qx,ABS(Qz))-Asp*S2*ZL(Qxp,ABS(Qzp))
	D1= -S2*sigmap*(1.5*VRTeTs(m)/SQRT(VRNeNs(m))*Qxp &
          /SQRT(1.+1.5*Beta*(Qxp2+Qzp2)))

!Decay normal:
	Funcx_RcdLLS= ZL(Qx,ABS(Qz))-Asp*S2*ZL(Qxp,ABS(Qzp)) -Aspp*S3*ZS(Qxdif,ABS(Qzdif))
	D1= -S2*sigmap*(1.5*VRTeTs(m)/SQRT(VRNeNs(m))*Qxp &
          /SQRT(1.+1.5*Beta*(Qxp2+Qzp2))) &
          + S1*S3*sigmapp*AA*SQRT(VRTeTs(m))*Qxdif &
          /SQRT(Qxdif**2+Qzdif**2+EpsMin) &
          /(1.+Beta*(Qxdif**2+Qzdif**2)/2.)**(1.5)

!LAproxS:
	Funcx_RcdLLS= ZL(Qx,ABS(Qz))-Asp*S2*ZL(Qxp,ABS(Qzp)) - Aspp*S3*ZSAprox(Qxdif,ABS(Qzdif))
	D1 = -S2*sigmap*(1.5*VRTeTs(m)/SQRT(VRNeNs(m))*Qxp &
          /SQRT(1.+1.5*Beta*(Qxp2+Qzp2))) &
          + S1*S3*sigmapp*AA*SQRT(VRTeTs(m))*Qxdif &
          /SQRT(Qxdif**2+Qzdif**2+EpsMin) &
          /(1.+Beta*(Qxdif**2+Qzdif**2)/2.)**(1.5)

!LAproxSD1:
	Funcx_RcdLLS= ZL(Qx,ABS(Qz))-Asp*S2*ZL(Qxp,ABS(Qzp)) - Aspp*S3*ZSAprox(Qxdif,ABS(Qzdif))
	D1 = -S2*sigmap*(1.5*VRTeTs(m)/SQRT(VRNeNs(m))*Qxp) &
          + S1*S3*sigmapp*AA*SQRT(VRTeTs(m))*Qxdif &
          /SQRT(Qxdif**2+Qzdif**2+EpsMin) &
          /(1.+Beta*(Qxdif**2+Qzdif**2)/2.)**(1.5)

!SaproxD1tambem:
	Funcx_RcdLLS= ZL(Qx,ABS(Qz))-Asp*S2*ZL(Qxp,ABS(Qzp)) - Aspp*S3*ZSAprox(Qxdif,ABS(Qzdif))
	D1AproxNoS= -S2*sigmap*(1.5*VRTeTs(m)/SQRT(VRNeNs(m))*Qxp &
          /SQRT(1.+1.5*Beta*(Qxp2+Qzp2))) &
          + S1*S3*sigmapp*AA*SQRT(VRTeTs(m))*Qxdif &
          /SQRT(Qxdif**2+Qzdif**2+EpsMin)

!SaproxD1tbMOD:
	Funcx_RcdLLS= ZL(Qx,ABS(Qz))-Asp*S2*ZL(Qxp,ABS(Qzp)) - Aspp*S3*ZSAprox(Qxdif,ABS(Qzdif))
	D1= -S2*sigmap*(1.5*VRTeTs(m)/SQRT(VRNeNs(m))*Qxp &
          /SQRT(1.+1.5*Beta*(Qxp2+Qzp2))) &
          + S1*S3*sigmapp*AA*SQRT(VRTeTs(m))*Qxdif &
          /SQRT(Qxdif**2+Qzdif**2 + 1) ! Modifiquei o + EpsMin para +1, porque pensei que podia estar ocorrendo algo ali

!
!
! Expressão para ZLAprox:
ZLAprox= SQRT(VRNeNs(m))*(1.+0.75*Q2*VRTeTs(m)/VRNeNs(m))

!
! Expressão para ZSAprox:
ZSAprox= Q*AA*SQRT(VRTeTs(m))