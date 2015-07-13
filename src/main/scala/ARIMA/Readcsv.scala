
package ARIMA

import scala.io._

object Readcsv {
  
  def read() = {
	val src = Source.fromFile("D://Rdata/test.csv")
	val iter = src.getLines().map(_.split(","))
	val data = iter.toList.map(x => x(2).toDouble).toArray
	data
  }
}