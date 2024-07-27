select ?v0 ?rating1 ?rating3 ?rating4 where
{
	?v0 <http://www4.wiwiss.fu-berlin.de/bizer/bsbm/v01/vocabulary/rating2> "6"^^<http://www.w3.org/2001/XMLSchema#integer> .
	?v0 <http://www4.wiwiss.fu-berlin.de/bizer/bsbm/v01/vocabulary/rating1> ?rating1.
	?v0 <http://www4.wiwiss.fu-berlin.de/bizer/bsbm/v01/vocabulary/rating3> ?rating3.
	?v0 <http://www4.wiwiss.fu-berlin.de/bizer/bsbm/v01/vocabulary/rating4> ?rating4.
}order by desc(?rating1+?rating3+?rating4) limit 5
