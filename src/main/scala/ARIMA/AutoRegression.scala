
/**
 * @author  Jian Wang
 * @version 1.0
 * @date    10/06/2015
 */

// package TimeSeriesForecast.ARIMA
package ARIMA

import breeze.linalg.{DenseVector, DenseMatrix, diag}
import stat.StatDenseVector

// Implementation of the Yule-Walker method by using Durbin-Levinson Algorithm to estimate the coefficient of AR model.

class AutoRegression (data : Array[Double], p : Int) extends Preprocessing(data : Array[Double]){
  
  private val phi = DenseMatrix.zeros[Double](maxOrder+1,maxOrder+1)
  private val v = DenseVector.zeros[Double](maxOrder+1)    
  
  def checkOrder : Boolean = p > maxOrder
  
  def durbinLevinson(data : Array[Double] = data, p : Int = p) : DenseVector[Double] = {
    v(0) = acvf(0)
	    
	for (t <- 1 to maxOrder) {
		var temp = 0.0
	    for (j <- 1 until t) {
	      temp += phi(t-1,j) * acvf(t-j)
	    }
	    phi(t,t) = (acvf(t) - temp) / v(t-1)
	    for (j <- 1 until t) phi(t,j) = phi(t-1,j) - phi(t,t) * phi(t-1,t-j)
	    v(t) = v(t-1)*(1. - phi(t,t)*phi(t,t))
	}
    
//	println(acvf)
//	println ("phi = " + phi)
	pacf = diag(phi) 
	pacf(0) = 1
	phi(p,1 to p).t 
  }
  
  def run() : AutoRegressionModel = {
    
//  a little bit different (maybe R's optimization?) with what in R with ar.yw() function. 
//  I strictly followed the book where sigma^2 = autocov(0) - (optimalphi dot (autocov(1),..autocov(p))))
// 	val sigsqr = acov(0) - (optphi dot acov(1 to p))		
//  println("sigmasqr = " + sigsqr) 
    
    if (checkOrder) sys.error("Order p is too large, larger than max order min(n-1, 20*log10(n)) or less than 1")
    
  	val optphi = durbinLevinson()
  	val sigma2 = v(p)
  	
//  	println("optphi =" + optphi.toArray.toList)
//  	println("sigma^2 = " + sigma2)  
  	
  	new AutoRegressionModel( y, mu, (optphi, sigma2))
  	
  }
}

class AutoRegressionModel (y : StatDenseVector, mu : Double , result : (DenseVector[Double], Double)) {
  val n = y.n
  val phi = result._1
  val sigma2 = result._2
  private val p = phi.length
  private def predictOne (newy : Array[Double], lastIDX : Int) : Double = {
    var sum = 0.0
    for (j <- 0 until p){
      sum += phi(j) * newy(lastIDX+1-p+j)      
    }
    sum
  }
  
  def predict (predictionLength: Int = 12) : Array[Double] = {
    val forecast = DenseVector.zeros[Double](predictionLength)
    val newy = new Array[Double](p + predictionLength)
    for (i <- 0 until p) newy(i) = y(n-p+i)
    var lastIDX = p-1
    for (i <- 0 until predictionLength){
      var next = 0.0
      next = predictOne(newy, lastIDX)
      forecast(i) = next
      newy(lastIDX+1) = next
      lastIDX += 1
    }  
    forecast.map(f => f + mu).toArray
  }
  
  def print_Result() : Unit = {
    println(result)
  }
}

object AutoRegression{
   def train(data: Array[Double], p: Int = 1): AutoRegressionModel = { 
     new AutoRegression(data, p).run()
  } 
}