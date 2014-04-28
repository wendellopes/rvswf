#-------------------------------------------------------------------------------
# Calculations must be performed for psi(n,x) and xi(n,x), where
# 1<=n<=N and x=x and x=Mx. Number of points: np=2*N*nx*ny*nz
#-------------------------------------------------------------------------------
# FRACAO CONTINUADA PELO METODO DE LENTZ
# f=bo+(a1/b1+)(a2/b2+)(a3/b3+)...(an/bn+...)
# input: vetores an e bn em que an=(a0,a1,...,aN) e bn=(b0,b1,...,bN)
# a0 may be zero or any number
#-------------------------------------------------------------------------------
# BASE NAMES
# lcf : Lentz Continued Fraction Evaluation
# ofc : Optical Force Calculations
# tst : Test Functions
# cmp : Comparison Functions
# aux : Auxiliary Functions
# FUNCTION NAMES
# cbld : Cylindrical Bessel Log Derivative
# cbri : Cylindrical Bessel Ratio Inverse (Barnett)
# cbrd : Cylindrical Bessel Ratio Direct
# sbld : Spherical Bessel Log Derivative
# sbri : Spherical Bessel Ratio Inverse (Barnett)
# sbrd : Spherical Bessel Ratio Direct
# rbld : Riccati-Bessel Log Derivative
# rbri : Riccati-Bessel Ratio Inverse (Barnett)
# rbrd : Riccati-Bessel Ratio Direct
# afsn : Auxiliary function Sn
