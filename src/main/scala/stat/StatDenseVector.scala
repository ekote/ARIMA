/**
 * @author  Jian Wang
 * @version 1.0
 * @date    10/06/2015
 */

package stat

import breeze.linalg._
import breeze.stats._
import breeze.numerics._
import breeze.stats.mean.reduce_Double
import breeze.stats.variance.reduceDouble

class StatDenseVector (X : Array[Double],unbiased: Boolean = false) extends DenseVector(X : Array[Double]){
  
	val n = X.length
	val mu = mean(X)
	val sum = X.sum
	val sig = variance(X)
	val den: Double = if (unbiased) n - 1.0 else n
	
	def this(X : DenseVector[Double]){
	  this(X.toArray)
	}
	
	// computer the covariance of two vector X and Y
	def cov(Y : DenseVector[Double]) : Double = {
	  val muY = mean(Y)
	  (this.dot(Y) -  sum * Y.sum / n) / den
	}
	
	// Compute the k-lag auto-covariance of this vector
	def acvf(k : Int = 1, demean : Boolean = true): Double = {
	  val nk = n - k
	  var sum = 0.0
	  for (i <- 0 until nk) {
	    if (demean) sum += (X(i) - mu) * (X(i+k) - mu)
	    else sum += X(i) * X(i+k)
	  }
	  sum / n
	}
			
	def acf(Type : String = "correlation", maxlag : Int = n - 1 ,demean : Boolean = true) : DenseVector[Double] = {
	  val acfv = DenseVector.zeros[Double](n)
	  for (t <- 0 until n){
	    if (demean) acfv(t) = acvf(t) else acfv(t) = acvf(t, demean = false)
	  } 
	  val cor = acfv / acfv(0)
	  if (Type == "correlation") cor(0 to maxlag)
	  else if (Type == "covariance") acfv(0 to maxlag)
	  else sys.error("Invalid Type. Please Type = correlation or covariance")
	}
	
	//Compute covariance ��(i,j)
	def k(i : Int, j : Int) : Double = acvf(i-j)
	
	// compute the Pearson Correlation of vector with another vector
	def corr (Y: StatDenseVector) : Double = cov (Y) / sqrt (sig * Y.sig )
	
	def toDenseVector : DenseVector[Double] = {
	  new DenseVector(X)
	}
}