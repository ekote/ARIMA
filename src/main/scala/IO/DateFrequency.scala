/**
 * @author  Jian Wang
 * @version 1.0
 * @date    10/06/2015
 */
package IO

import Readcsv._
import breeze.linalg.{DenseVector}

// this object is to change the daily data to the specific data, weekly if set steps to 7 and monthly if set steps to 30
object DateFrequency {
  
  def FromDaily(data: Array[Double],steps: Int = 7): Array[Double] = {
    val n = data.length
    val y = Array.fill[Double](n/steps)(0.0)
    
    var sum = 0.0 
    for (i <- 0 until n) {
      sum += data(i)
      if ((i+1) % steps == 0) {
        y((i+1)/steps-1) = sum
        sum = 0.0
      }
    }
    y
  }
  
  
}