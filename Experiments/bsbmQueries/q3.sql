select ?v0 ?v1 ?v2 ?v3 ?v4  where
{
	?v0 <http://purl.org/dc/elements/1.1/title> ?v1 .
	?v0 <http://purl.org/dc/elements/1.1/publisher> <http://www4.wiwiss.fu-berlin.de/bizer/bsbm/v01/instances/dataFromRatingSite11/RatingSite11> .
	?v0 <http://purl.org/stuff/rev#text> ?v4 .
	?v2 <http://purl.org/stuff/rev#reviewer> ?v3 .
	?v2 <http://purl.org/dc/elements/1.1/date> "2008-04-16"^^<http://www.w3.org/2001/XMLSchema#date> .
	?v2 <http://purl.org/stuff/rev#text> ?v4.
	?v0 <http://www4.wiwiss.fu-berlin.de/bizer/bsbm/v01/vocabulary/rating1> ?rank1 .
	?v2 <http://www4.wiwiss.fu-berlin.de/bizer/bsbm/v01/vocabulary/rating1> ?rank2 .
}order by (-?rank1-?rank2) limit 5
