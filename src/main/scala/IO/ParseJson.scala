package IO

import scala.io.{ Source }
import org.json4s._
import org.json4s.native.JsonMethods._

object ParseJson {
  def parseJson(path: String, startdate: String, enddate: String) = {
    implicit val formats = DefaultFormats
    val src = Source.fromFile(path)
    val json = try src.mkString finally src.close()
    val jobject = parse(json, useBigDecimalForDouble = true)

    val rows = for {
      JObject(queries) <- jobject
      JField("rows", JArray(rows)) <- queries
    } yield rows

    val data = rows(0).map(v => v.extract[List[String]]).map(v => (v(0), v(1).toDouble)).toArray
    data.filter(v => (v._1 >= startdate) && (v._1 < enddate))
  }

}