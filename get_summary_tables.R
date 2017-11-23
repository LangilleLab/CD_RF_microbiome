# Loop through RF summary objects (.rda files) and make summary table for disease and treatment response

original_colnames <-  c( "file" , "num.features" , "mtry" , "num.trees" , "LOOCV.Accuracy" , "Model.OOB.error" , "Median.Rand.OOB.error" , "P.value" )

summary_table <- data.frame( matrix( c( NA ) , nrow = 0 , ncol = 11 ) )

colnames(summary_table) <- c( "trait" , "sequencing" , "category"  , original_colnames )

sequencing <- c( "16S" , "MGS" )

# Note that KOs are not included since they are run with a lower number of trees below
categories <- c( "phylum" , "class" , "order" , "family" , "genus" , "species", "otus" , "strains" , "modules" , "pathways" , "KOs" )
traits <- c( "disease" , "response" )

for ( s in sequencing ) {

  for ( c in categories ) {

    for ( t in traits ) {

      if ( c == "KOs") {

        rda <- paste( "RF_obj_out/KOs/" , s  , sep ="" )

      } else {

        rda <- paste( "RF_obj_out/" , s  , sep ="" )
      }

      rda <- paste( rda , c , t , "RF_out.rda" , sep="_" ) 

      if ( file.exists( rda ) ) {

        attach( rda )

        reordered_summary <- rf_out$summary

        reordered_summary$sequencing <- s
        reordered_summary$category <- c
        reordered_summary$trait <- t

        reordered_summary <- reordered_summary[ , c( "trait" , "sequencing" , "category"  , original_colnames ) ]

        summary_table <- rbind( summary_table , reordered_summary )

        detach()

      }

    }

  }

}

write.table( x=summary_table , file="RF_summary_table.txt" , row.names=FALSE , quote = FALSE , sep="\t" )

