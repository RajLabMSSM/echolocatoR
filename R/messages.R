


# Keep separate bc the emojis can mess up the scripts....
messages <- function(name){
  message_list <- list(locus="🦇 🦇 🦇 🦇 🦇 🦇",
                       query="\n------------------ Step 1: Query 🔍 ---------------",
                       LD = "\n--- Step 2: Extract Linkage Disequilibrium 🔺 --",
                       filter = "\n-------------- Step 3: Filter SNPs 🚰 -------------",
                       finemap="\n-------- Step 4: Fine-map 🔊 --------",
                       visualize="\n--------------- Step 5: Visualize 📊--------------"
                       )
  message(message_list[[name]])
}
