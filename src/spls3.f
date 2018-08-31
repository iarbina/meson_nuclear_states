      SUBROUTINE SPLS3(X,Y,N,XI,FI,M,Q,AU,IGO,ISPL)                     
!                                                                       
!     ******************************************************************
!                                                                       
!     CUBIC SPLINE, STARTING WITH ZERO SECOND DERIVATIVES AT THE        
!     BOUNDARIES OF THE APPROXIlopez1991!IMATION INTERVAL.                         
!                                                                       
!     IGO = 0      BUILD UP SPLINE ONLY.                                
!     IGO = 1      BUILD UP SPLINE AND INTERPOLATE FOR M POINTS.        
!     IGO = 2      BUILD UP SPLINE AND COMPUTE DERIVATIVES AT M POINTS. 
!                                                                       
!     REAL*8 VERSION.        J.GALONSKA, 15.12.1971           
!     ******************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)                                       
      DIMENSION X(1),Y(1),XI(1),Q(1),AU(1),FI(1)
!                                                                       
      ZERO=0.D0                                                         
      THREE=3.D0                                                        
      SIX=6.D0                                                          
      FACT=0.1666666666667D0                                            
      IF (ISPL.NE.0)  GO TO 30                                          
!                                                                       
!                                                                       
      AU(1) = ZERO                                                      
      AU(N) = ZERO                                                      
      Q(1) = ZERO                                                       
      HK = X(2) - X(1)                                                  
      YSAVE = (Y(2)-Y(1)) / HK                                          
      AUX = ZERO                                                        
      NN = N - 1                                                        
      DO 10  K = 2,NN                                                   
      HX = X(K+1) - X(K-1)                                              
      DIVQ = (HK*Q(K-1)+HX+HX)                                          
      HK = X(K+1) - X(K)                                                
      YK = (Y(K+1)-Y(K)) / HK                                           
      Q(K) = - HK / DIVQ                                                
      AU(K) = (SIX*(YK-YSAVE)-AUX) / DIVQ                               
      YSAVE = YK                                                        
      AUX = AU(K) * HK                                                  
   10 CONTINUE                                                          
!                                                                       
      NN2 = NN + 2                                                      
      DO 20  KK = 2,NN                                                  
      K = NN2 - KK                                                      
   20 AU(K) = Q(K) * AU(K+1) + AU(K)                                    
!                                                                       
      IF (IGO.EQ.0)  RETURN                                             
!                                                                       
!     ******************************************************************
!                                                                       
!     INTERPOLATION OR COMPUTATION OF DERIVATIVES.                      
!                                                                       
!     IGO = 1      INTERPOLATE FOR M POINTS.                            
!     IGO = 2      COMPUTE DERIVATIVES AT M POINTS.                     
!                                                                       
!     ******************************************************************
!                                                                       
   30 DO 100  J = 1,M                                                   
      IF (X(1).GT.XI(J))  THEN                                         
      M1=1
      M2=2
      GO TO 90
      ELSE
      ENDIF
      IF (XI(J).GT.X(N))  THEN                                       
      M1=N-1
      M2=N
      GO TO 90
      ELSE
      ENDIF
      M1 = 1                                                            
      M2 = N                                                            
   50 M3 = (M2+M1)/2                                                    
      IF (XI(J).GE.X(M3))  GO TO 70                                     
      M2 = M3                                                           
      GO TO 80                                                          
   70 M1 = M3                                                           
   80 IF (M1+1-M2.NE.0)  GO TO 50                                       
   90 DIJ = X(M2) - XI(J)                                               
      DIM1J = X(M1) - XI(J)                                             
      HI = X(M2) - X(M1)                                                
      HI2 = HI * HI                                                     
      IF (IGO.GE.2)  GO TO 95                                           
      DIJ3 = DIJ * DIJ * DIJ                                            
      DIM1J3 = DIM1J * DIM1J * DIM1J                                    
      FI(J) = FACT * (AU(M1)*DIJ3-AU(M2)*DIM1J3+(SIX*Y(M1)-HI2*AU(M1))
     &        *DIJ-(SIX*Y(M2)-HI2*AU(M2))*DIM1J) / HI                   
      GO TO 100                                                         
   95 FI(J) = FACT * (THREE*(AU(M2)*DIM1J*DIM1J-AU(M1)*DIJ*DIJ)
     &       -SIX*(Y(M1)-Y(M2))+HI2*(AU(M1)-AU(M2))) / HI  
  100 CONTINUE               
      RETURN                                                          
!                                                                       
! 110 M1 = 1                                                            
!     M2 = 2                                                            
!     GO TO 90                                                          
!                                                                       
! 120 M1 = N - 1                                                        
!     M2 = N                                                            
!     GO TO 90                                                          
!                                                                       
      END SUBROUTINE SPLS3
