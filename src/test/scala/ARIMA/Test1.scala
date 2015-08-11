package ARIMA

import org.specs2._
import IO.Readcsv._

class Test1 extends Specification {
  def is = s2"""

 This is a specification to check ARIMA model based on R function arma and forecast in package itsmr
  'Note: itsmr calls arima in R stats'
 Test R code is as following:
  
  library(itsmr)
  data <- read.csv("D:/workspace/ARIMA/src/test/resources/timeseries_ppi.csv")
  attach(data)
  y <- ppi
  e = Resid(y,xv = c("diff", 1))
  a = autofit(e)
  # a = arma(e,3,2)
  # which is equvalent to call arima(e,order = c(3,0,2)) in R stats
  forecast(y,xv = c("diff",1), a, h = 12)
  
 The length of phi should be 3 $phiSize
 The length of theta should be 2 $thetaSize
 phi1 is around 0.1813941 +- 0.1 $phi1
 phi2 is around -0.6311208 +- 0.1 $phi2
 phi3 is around 0.6235352 +- 0.1 $phi3
 theta1 is around 0.3484732 +- 0.1 $theta1
 theta2 is around 0.8860058 +- 0.1 $theta2
 aicc is around 380.704 +- 0.1 $aicc
 
 Now test forecasting
 predictions, Standard errors (not included if there is a log transform)
 LowerBound, UpperBound of the predictions
 
 $pred1
 $pred2
 $pred3
 $pred4
 $pred5
 $pred6
 $pred7
 $pred8
 $pred9
 $pred10
 $pred11
 $pred12
 
 $sqrtMSE1
 $sqrtMSE2
 $sqrtMSE3
 $sqrtMSE4
 $sqrtMSE5
 $sqrtMSE6
 $sqrtMSE7
 $sqrtMSE8
 $sqrtMSE9
 $sqrtMSE10
 $sqrtMSE11
 $sqrtMSE12
 
 $l1
 $l2
 $l3
 $l4
 $l5
 $l6
 $l7
 $l8
 $l9
 $l10
 $l11
 $l12
 
 $r1
 $r2
 $r3
 $r4
 $r5
 $r6
 $r7
 $r8
 $r9
 $r10
 $r11
 $r12
 
 
"""

  val x = csvReader("D://Data/sharefolder/testdata.csv")
  val model = ARIMA.train(x, order = (3, 1, 2), xv = List("diff", "1"), method = "MLE", demean = true)
  //  val model = AutoFit.autofit(x)
  def phiSize = model.phiHat.toArray must have length (3)
  def thetaSize = model.thetaHat.toArray must have length (2)
  def phi1 = model.phiHat.toArray.apply(0) must beCloseTo(0.1813941, 0.1)
  def phi2 = model.phiHat.toArray.apply(1) must beCloseTo(-0.6311208, 0.1)
  def phi3 = model.phiHat.toArray.apply(2) must beCloseTo(0.6235352, 0.1)
  def theta1 = model.thetaHat.toArray.apply(0) must beCloseTo(0.3484732, 0.1)
  def theta2 = model.thetaHat.toArray.apply(1) must beCloseTo(0.8860058, 0.1)
  def aicc = model.aicc must beCloseTo(380.704, 0.1)

  val predictions = model.predict(12)
  val pred = predictions._1
  val sqrtMSE = predictions._2
  val lowerBound = predictions._3
  val upperBound = predictions._4

  def pred1 = pred(0) must beCloseTo(102.6415, 0.1)
  def pred2 = pred(1) must beCloseTo(102.3812, 0.1)
  def pred3 = pred(2) must beCloseTo(103.0154, 0.1)
  def pred4 = pred(3) must beCloseTo(103.2054, 0.1)
  def pred5 = pred(4) must beCloseTo(103.0609, 0.1)
  def pred6 = pred(5) must beCloseTo(103.6938, 0.1)
  def pred7 = pred(6) must beCloseTo(104.4019, 0.1)
  def pred8 = pred(7) must beCloseTo(104.4244, 0.1)
  def pred9 = pred(8) must beCloseTo(104.7598, 0.1)
  def pred10 = pred(9) must beCloseTo(105.6315, 0.1)
  def pred11 = pred(10) must beCloseTo(105.9756, 0.1)
  def pred12 = pred(11) must beCloseTo(106.0805, 0.1)

  def sqrtMSE1 = sqrtMSE(0) must beCloseTo(0.7205607, 0.1)
  def sqrtMSE2 = sqrtMSE(1) must beCloseTo(1.31697, 0.1)
  def sqrtMSE3 = sqrtMSE(2) must beCloseTo(1.88976, 0.1)
  def sqrtMSE4 = sqrtMSE(3) must beCloseTo(2.482266, 0.1)
  def sqrtMSE5 = sqrtMSE(4) must beCloseTo(3.027963, 0.1)
  def sqrtMSE6 = sqrtMSE(5) must beCloseTo(3.49925, 0.1)
  def sqrtMSE7 = sqrtMSE(6) must beCloseTo(3.952317, 0.1)
  def sqrtMSE8 = sqrtMSE(7) must beCloseTo(4.392986, 0.1)
  def sqrtMSE9 = sqrtMSE(8) must beCloseTo(4.783067, 0.1)
  def sqrtMSE10 = sqrtMSE(9) must beCloseTo(5.142448, 0.1)
  def sqrtMSE11 = sqrtMSE(10) must beCloseTo(5.501305, 0.1)
  def sqrtMSE12 = sqrtMSE(11) must beCloseTo(5.837452, 0.1)

  def l1 = lowerBound(0) must beCloseTo(101.2292, 0.1)
  def l2 = lowerBound(1) must beCloseTo(99.79994, 0.1)
  def l3 = lowerBound(2) must beCloseTo(99.31149, 0.1)
  def l4 = lowerBound(3) must beCloseTo(98.34019, 0.1)
  def l5 = lowerBound(4) must beCloseTo(97.12608, 0.1)
  def l6 = lowerBound(5) must beCloseTo(96.83526, 0.1)
  def l7 = lowerBound(6) must beCloseTo(96.65535, 0.1)
  def l8 = lowerBound(7) must beCloseTo(95.81411, 0.1)
  def l9 = lowerBound(8) must beCloseTo(95.38495, 0.1)
  def l10 = lowerBound(9) must beCloseTo(95.55234, 0.1)
  def l11 = lowerBound(10) must beCloseTo(95.19303, 0.1)
  def l12 = lowerBound(11) must beCloseTo(94.63912, 0.1)

  def r1 = upperBound(0) must beCloseTo(104.0538, 0.1)
  def r2 = upperBound(1) must beCloseTo(104.9625, 0.1)
  def r3 = upperBound(2) must beCloseTo(106.7193, 0.1)
  def r4 = upperBound(3) must beCloseTo(108.0707, 0.1)
  def r5 = upperBound(4) must beCloseTo(108.9957, 0.1)
  def r6 = upperBound(5) must beCloseTo(110.5523, 0.1)
  def r7 = upperBound(6) must beCloseTo(112.1484, 0.1)
  def r8 = upperBound(7) must beCloseTo(113.0346, 0.1)
  def r9 = upperBound(8) must beCloseTo(114.1346, 0.1)
  def r10 = upperBound(9) must beCloseTo(115.7107, 0.1)
  def r11 = upperBound(10) must beCloseTo(116.7581, 0.1)
  def r12 = upperBound(11) must beCloseTo(117.5219, 0.1)

}