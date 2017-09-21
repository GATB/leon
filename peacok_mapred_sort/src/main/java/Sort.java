import java.io.IOException;
import java.lang.Runtime;
import java.lang.Process;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.OutputStream;
import static java.lang.Math.toIntExact;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.conf.Configurable;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.mapreduce.Mapper;
import org.apache.hadoop.mapreduce.Reducer;
import org.apache.hadoop.mapreduce.Partitioner;
import org.apache.hadoop.mapreduce.lib.partition.TotalOrderPartitioner;
import org.apache.hadoop.mapreduce.lib.partition.InputSampler;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.mapreduce.lib.input.TextInputFormat;
import org.apache.hadoop.mapreduce.lib.output.TextOutputFormat;
//import org.apache.hadoop.util;
import java.util.Map;
import java.util.AbstractMap;
import java.util.TreeMap;
import java.net.URI;


/*class MutableInt {
	  long value = 0; 
	  
	  MutableInt(long val){
		  value = val;
	  }
	  
	  public void increment () { ++value;      }
	  public long  get ()       { return value; }
}
*/

class InfosReducers{
	
	int numActualReducerID;
	//int nbReducersNeeded;
	long nbTasksLeftForActualReducer;
	
	InfosReducers(int numActualReducerID, long nbTasksLeftForActualReducer){
		
		this.numActualReducerID = numActualReducerID;
		this.nbTasksLeftForActualReducer = nbTasksLeftForActualReducer;
	}
	
	public int getNumActualReducerID() {
		return numActualReducerID;
	}
	
	/*public int getNbReducersNeeded() {
		return nbReducersNeeded;
	}*/
	
	public long getNbTasksLeftForActualReducer() {
		return nbTasksLeftForActualReducer;
	}
	
	//give the id of next reducer
	//actualize nb of tasks remaining to attribute to actual reducer
	//if no place left, next tasks will be attributed to the followinf reducer
	public int getNextReducerID(long limitReducerCapacity) {
		
		if (nbTasksLeftForActualReducer > 0) {
			--nbTasksLeftForActualReducer;
		}
		else {
			++numActualReducerID;
			//--nbReducersNeeded;
			nbTasksLeftForActualReducer = limitReducerCapacity;
			
		}
		return numActualReducerID;
	}
}

public class Sort {
	
	//static Map<String, MutableInt> histoKeys = new TreeMap<String, MutableInt>();
	
	static int numReducers = 0;
	static Map<String, Long> histoTreeMap = new TreeMap<String, Long>();
	static Map<String, Integer> reducerTreeMap = new TreeMap<String, Integer>();
	static Map<String, InfosReducers> infosReducersTreeMap = new TreeMap<String, InfosReducers>();
	
	
	public static class HistoMapper
    extends Mapper<Object, Text, Text, IntWritable>{
	  
		public void map(Object key, Text value, Context context
			  ) throws IOException, InterruptedException {
		  
		  String tokens[] = value.toString().split(";");

		  Text anchor = new Text(tokens[0]);
		  Text readInfos = new Text(tokens[1]);
		  
		  //counting occurence of each keys for further partitionning
		  //String anchorString = anchor.toString();
		  //MutableInt count = histoKeys.get(anchorString);
		  //if (count == null) {
		//	  histoKeys.put(anchorString, new MutableInt(1));
		  //}
		  //else {
		  //    count.increment();
		  //}
		  
		  //++totalNbElements;
		  
		  
		  
		  context.write(anchor, new IntWritable(1));
	  }
	}
	
	public static class HistoReducer
      extends Reducer<Text,IntWritable,Text,LongWritable> {
   
		
	  public void reduce(Text key, Iterable<IntWritable> values,
                      Context context
                      ) throws IOException, InterruptedException {
		  
		  
		  long sum = 0;
		  
		  for (IntWritable val : values){
			  sum += val.get();
		  }
		  
		  //String readsInfosString = "\n";
		  //for (Text val : values){
			  
			//  readsInfosString += val.toString();
			//  readsInfosString += "\n";
		  //}
		  
		  //Text readsInfos = new Text(readsInfosString);
		  
		  
		  context.write(key, new LongWritable(sum));
	  }
	}
	
	public static class SortMapper
       extends Mapper<Object, Text, Text, Text>{
	  
		public void map(Object key, Text value, Context context
			  ) throws IOException, InterruptedException {
		  
		  String tokens[] = value.toString().split(";");

		  Text anchor = new Text(tokens[0]);
		  Text readInfos = new Text(tokens[1]);
		  
		  //counting occurence of each keys for further partitionning
		 /* String anchorString = anchor.toString();
		  MutableInt count = histoKeys.get(anchorString);
		  if (count == null) {
			  histoKeys.put(anchorString, new MutableInt());
		  }
		  else {
		      count.increment();
		  }
		  
		  ++totalNbElements;
		  */
		  
		  context.write(anchor, readInfos);
	  }
  }
	
	
	public static class SortPartitioner
		extends Partitioner<Text, Text>
		implements Configurable
	{
		private Configuration configuration;
		static Map<String, Integer> reducerTreeMap =  new TreeMap<String, Integer>();
		static Map<String, InfosReducers> infosReducersTreeMap = new TreeMap<String, InfosReducers>();
		//boolean reducerTreeMapSet = false;
		//Map<String, Integer> keyReducer = new TreeMap<String, Integer>();
		int partitionVersion;
		long limitReducerCapacity;
		
		/*public SortPartitioner()
		{			
			//create the map of key/reducer 
			//System.out.println(Sort.numReducers);
			//float boundPortion = totalNbElements/Sort.numReducers+1;
			/*long currentCount = 0;
			int reducerID = 0;
			
			for(Map.Entry<String,MutableInt> entry : histoKeys.entrySet()) {
				  String key = entry.getKey().toString();
				  long value = entry.getValue().get();

				   
				  currentCount += value;
				  
				  
				  if (currentCount <= boundPortion)
				  {
					  keyReducer.put(key, reducerID);
				  }
				  else
				  {

					++reducerID;
					currentCount = 0;
					  
					keyReducer.put(key, reducerID);
				  }
				}
			*/
			//super();
			//String sreducerTreeMap = configuration.get("sreducerTreeMap");
			
			//String reducerTreeMapKeysValues[] = sreducerTreeMap.split(";");
			
			// for (String reducerTreeMapKeyValue : reducerTreeMapKeysValues){
				 
			//	 String tokens[] = reducerTreeMapKeyValue.split(":");
				 
			//	 reducerTreeMap.put(tokens[0], Integer.parseInt(tokens[1]));
			 // }
		/*	  
			  String treeMapFilePath = configuration.get("treeMapFilePath");
			
			try{
		        Path pt=new Path(treeMapFilePath);
		        FileSystem fs = FileSystem.get(configuration);
		        BufferedReader br=new BufferedReader(new InputStreamReader(fs.open(pt)));
		        String line;
		        line=br.readLine();
		        while (line != null){
		                System.out.println(line);
		                
		                String tokens[] = line.split(":");
		                
		                String anchor = tokens[0];
		                int reducer = Integer.parseInt(tokens[1]);
		                
		                reducerTreeMap.put(anchor, reducer);
		                line=br.readLine();
		                
		        }
			}catch(Exception e){
				System.out.println("\t\tERROR WHILE READING");
			}
		}*/
		
		private void setReducerTreeMap(){
			
			String treeMapFilePath = configuration.get("treeMapFilePath");
			
			try{
		        Path pt=new Path(treeMapFilePath);
		        FileSystem fs = FileSystem.get(configuration);
		        BufferedReader br=new BufferedReader(new InputStreamReader(fs.open(pt)));
		        String line;
		        line=br.readLine();
		        while (line != null){
		                //System.out.println(line);
		                
		                String tokens[] = line.split(":");
		                
		                String anchor = tokens[0];
		                int reducer = Integer.parseInt(tokens[1]);
		                
		                reducerTreeMap.put(anchor, reducer);
		                line=br.readLine();
		                
		        }
			}catch(Exception e){
				System.out.println("\t\tERROR WHILE READING");
			}
			
			
			/*String reducerTreeMapKeysValues[] = sreducerTreeMap.split(";");
	
			 for (String reducerTreeMapKeyValue : reducerTreeMapKeysValues){
				 String tokens[] = reducerTreeMapKeyValue.split(":");
				 reducerTreeMap.put(tokens[0], Integer.parseInt(tokens[1]));
			 }*/
		}
		
		
		@Override
	     public int getPartition(Text key, Text value, int numReduceTasks)
	      {
			
			/*if (!reducerTreeMapSet){
				setReducerTreeMap();
			}*/
			
	         String anchorString = key.toString();
	         
	         if(numReduceTasks == 0)
	         {
	            return 0;
	         }
	         
	         if(partitionVersion < 4) {
	        	 return reducerTreeMap.get(anchorString);
	         }
	         else {
	        	 InfosReducers ir = infosReducersTreeMap.get(anchorString);
	        	 return ir.getNextReducerID(limitReducerCapacity);
	         }
	         
	         //return 1;
	      }
		
		@Override
	      public void setConf(Configuration configuration) {
	         this.configuration = configuration;
	         
	         String treeMapFilePath = configuration.get("treeMapFilePath");
	         int partitionVersion = Integer.parseInt(configuration.get("partitionVersion"));
	         long limitReducerCapacity = Long.parseLong(configuration.get("limitReducerCapacity"));
				
	         try{
	        	 
		         Path pt=new Path(treeMapFilePath);
			     FileSystem fs = FileSystem.get(configuration);
			     BufferedReader br=new BufferedReader(new InputStreamReader(fs.open(pt)));
			     String line;
			     line=br.readLine();
		        
		     	if (partitionVersion < 4) {
	        	 
				
			        while (line != null){
			                //System.out.println(line);
			                
			                String tokens[] = line.split(":");
			                
			                String anchor = tokens[0];
			                int reducer = Integer.parseInt(tokens[1]);
			                
			                reducerTreeMap.put(anchor, reducer);
			                line=br.readLine();
			                
			        }
				}
		     	else {

		     		while (line != null){
		     			//System.out.println(line);
				                
		     			String tokens[] = line.split(":");
	        			 
		     			String anchor = tokens[0];
		     			int numFirstReducer = Integer.parseInt(tokens[1]);
		     			long nbTasksLeftForActualReducer = Long.parseLong(tokens[2]);
	        				
		     			InfosReducers ir = new InfosReducers(numFirstReducer, nbTasksLeftForActualReducer);

		     			infosReducersTreeMap.put(anchor, ir);
		     			line=br.readLine();
				                
		     		}
		     	}
	         }catch(Exception e){
        		 System.out.println("\t\tERROR WHILE READING"); 
	         }
	         
	      }
	 
	      @Override
	      public Configuration getConf() {
	         return configuration;
	      }
		
	}
	
	
  public static class SortReducer
       extends Reducer<Text,Text,Text,Text> {
    
	  public void reduce(Text key, Iterable<Text> values,
                       Context context
                       ) throws IOException, InterruptedException {
		  
		  Text anchor = key;
		  
		  
		  String readsInfosString = "\n";
		  for (Text val : values){
			  
			  readsInfosString += val.toString();
			  readsInfosString += "\n";
		  }
		  
		  Text readsInfos = new Text(readsInfosString);
		  
		  context.write(anchor, readsInfos);
    }
  }
  
  /*
   *version 1 
   *try to balance charge on remaining reducers
   *(the first one will probably be the more charged reducer)
   */
  public static void fillReducerTreeMapV1(long totalNbElements, int numReduceTasks) {
	  
	  long boundPortion = totalNbElements/numReducers+1;
	  long nbRemainingElements = totalNbElements;
	  int nbRemainingReducers = numReducers;
	  long limit = boundPortion;
	  long distance = boundPortion;
	  long newDistance = boundPortion;
	  long currentCount = 0;
	  int reducerID = 0;
		
	  for(Map.Entry<String,Long> entry : histoTreeMap.entrySet()) {
		  String anchor = entry.getKey().toString();
		  long occurence = entry.getValue();
	  
		  currentCount += occurence;
    	
		  reducerTreeMap.put(anchor, reducerID);
				  
		  if (currentCount > limit){

			  ++reducerID;
				
			  nbRemainingElements -= currentCount;
			  nbRemainingReducers -= 1;
			  boundPortion = nbRemainingElements/nbRemainingReducers+1;
			  currentCount = 0;
			  limit = boundPortion;
				
		  }
	  }
  }
  
  /*
   * version 2
   * balance by portions :
   * the first one will at least work on the first n ones
   * the second from n to 2n... etc
   * if one portion reaches the limit of the next reducer,
   * then the next one won't work at all...
   */
  public static void fillReducerTreeMapV2(long totalNbElements, int numReduceTasks) {
	  
	  long boundPortion = totalNbElements/numReducers+1;
	  long nbRemainingElements = totalNbElements;
	  int nbRemainingReducers = numReducers;
	  long limit = boundPortion;
	  long distance = boundPortion;
	  long newDistance = boundPortion;
	  long currentCount = 0;
	  int reducerID = 0;
		
	  for(Map.Entry<String,Long> entry : histoTreeMap.entrySet()) {
		  
		  String anchor = entry.getKey().toString();
		  long occurence = entry.getValue();
	  
	      currentCount += occurence;	    	
	      reducerTreeMap.put(anchor, reducerID);
				
	      if (currentCount > limit){
	
	    	  ++reducerID;
	    	  limit += boundPortion;
					
	      }
	  }
  }
  
  /*
   * version 3
   * take the next set in the reducer if the distance to
   * bound portion get smaller, and then reajust the perfect
   * bound
   * else stop
   */
  public static void fillReducerTreeMapV3(long totalNbElements, int numReduceTasks) {
	  
	  long boundPortion = totalNbElements/numReducers+1;
	  long nbRemainingElements = totalNbElements;
	  int nbRemainingReducers = numReducers;
	  long limit = boundPortion;
	  long distance = boundPortion;
	  long newDistance = boundPortion;
	  long currentCount = 0;
	  int reducerID = 0;
		
	  for(Map.Entry<String,Long> entry : histoTreeMap.entrySet()) {
		  
		  String anchor = entry.getKey().toString();
		  long occurence = entry.getValue();
	  
	      currentCount += occurence;
	
	      newDistance = Math.abs(boundPortion - currentCount);
	
				     
	      if (newDistance <= distance)
	      { 
	    	  reducerTreeMap.put(anchor, reducerID);
	    	  distance = newDistance;
	      }
	      else{
	    	  if (reducerID < numReduceTasks-1){
	    		  ++reducerID;
	    	  }
	    	  reducerTreeMap.put(anchor, reducerID);
	    	  nbRemainingElements -= currentCount;
	    	  nbRemainingReducers -= 1;
	    	  if (nbRemainingReducers > 0){
	    		  boundPortion = nbRemainingElements/nbRemainingReducers+1;
	    		  System.out.println("test - new boundPortion : " + boundPortion);
	    	  }
						
	    	  currentCount = 0;
	    	  distance = boundPortion;
	    	  newDistance = 0;
	      }
	  }
  }
  
  //V4 :
 // cut exactly charges between reducers
  //several red may have the same key
  //will have to remove first line of reducers resulting file, receiving
  // the same key that the previous one
public static void fillInfosReducersTreeMapV4(long totalNbElements, int numReduceTasks) {
	  
	  long boundPortion = totalNbElements/numReducers+1;
	  long nbRemainingElements = totalNbElements;
	  int nbRemainingReducers = numReducers;
	  long limit = boundPortion;
	  long distance = boundPortion;
	  long newDistance = boundPortion;
	  long currentCount = 0;
	  int reducerID = 0;
		
	  for(Map.Entry<String,Long> entry : histoTreeMap.entrySet()) {
		  String anchor = entry.getKey().toString();
		  long occurence = entry.getValue();
		  
		  int numFirstReducer;
		  long nbTasksLeftForActualReducer;
			
		  currentCount += occurence;
		  
		  numFirstReducer = reducerID;
				  
		  if (currentCount > limit){
			  
			 long remainingCount = currentCount - limit;
			 //number of red needed are 1 + the number of red for the exceeded data (remainingCount/limit)
			// int nbReducersNeeded = 1 + ((int) (remainingCount/limit));
			 int nbReducersNeeded = 1 + (int) Math.ceil((double) remainingCount / limit);
			 nbTasksLeftForActualReducer = occurence - remainingCount;
				
			 currentCount = remainingCount;
			 
			 
			//we increment reducerID wirh the nb red needed minus 1
			 //if the last one isn't full
			 if ((remainingCount % limit) != 0)
			 {
				 reducerID += nbReducersNeeded -1;
			 }
			 else
			 {
				 reducerID += nbReducersNeeded;
			 }
			 
			 
			 //old version, adjusting limit. Not necessary anymore (normaly)
			  //nbRemainingElements -= currentCount;
			  //nbRemainingReducers -= 1;
			  //boundPortion = nbRemainingElements/nbRemainingReducers+1;
			  //currentCount = 0;
			  //limit = boundPortion;
				
		  }
		  
		  if (currentCount == limit){
			  
			  //numFirstReducer = reducerID;
			  nbTasksLeftForActualReducer = occurence;
			  
			  ++reducerID;
			  currentCount = 0;
		  }
		  
		  else {
			  //numFirstReducer = reducerID;
			  nbTasksLeftForActualReducer = occurence;
		  }
		  
		  InfosReducers ir = new InfosReducers(numFirstReducer, nbTasksLeftForActualReducer);
		  infosReducersTreeMap.put(anchor, ir);
	  }
  }
  
  public static void main(String[] args) throws Exception {
	  
	if (args.length != 4){
		System.out.println("cmd inputPath outputPath nbElements nbReducers");
		System.exit(1);;
	}
	
	int partitionVersion = 1;
	  
	long totalNbElements = Long.parseLong(args[2]);  
    Configuration conf = new Configuration();
    int numReduceTasks = Integer.parseInt(args[3]);
    
    String arg0PartPath[] = args[0].split("/");
	String hdfsBasePath = arg0PartPath[0]+"/"+arg0PartPath[1]+"/"+arg0PartPath[2];
	//System.out.println("hdfsBasePath : "+hdfsBasePath);
	FileSystem hdfs = FileSystem.get( new URI(hdfsBasePath), conf);
    
    //create input path string
    String inputFilePath = "";
    String inputPath = "";
    if (null != inputFilePath && inputFilePath.length() > 0 )
    {
        int endIndex = inputFilePath.lastIndexOf("/");
        if (endIndex != -1)  
        {
        	inputPath = inputFilePath.substring(0, endIndex); // not forgot to put check if(endIndex != -1)
        }
    }
    
    
    //HISTOMAP
     
    String histosFilesPath = args[1] + ".tmp";
    
    //Create Job
    Job job1 = Job.getInstance(conf, "Histo");  
    job1.setJarByClass(Sort.class);
    
    //Number of Reducer tasks.
    job1.setNumReduceTasks(numReduceTasks);
    numReducers = job1.getNumReduceTasks();
    //System.out.println(numReducers);
    
    // File Input and Output paths
    FileInputFormat.addInputPath(job1, new Path(args[0]));
    FileOutputFormat.setOutputPath(job1, new Path(histosFilesPath));
    
    //Set Mapper class and Output format for key-value pair
    job1.setMapperClass(HistoMapper.class);
    job1.setMapOutputKeyClass(Text.class);
    job1.setMapOutputValueClass(IntWritable.class);
    
    
    //set partitioner statement
    //job1.setPartitionerClass(SortPartitioner.class);
    
    
    //Set Reducer class and Input/Output format for key-value pair
    job1.setReducerClass(HistoReducer.class);
    
    //Input and Output format for data
    job1.setInputFormatClass(TextInputFormat.class);
    job1.setOutputFormatClass(TextOutputFormat.class);
    job1.setOutputKeyClass(Text.class);
    job1.setOutputValueClass(LongWritable.class);
    
    job1.waitForCompletion(true);
    
    // CONSTRUCT HISTOGRAM
    
  //create the map of key/reducer 
    
    //merge reducers files
    String command = "hdfs dfs -getmerge " + histosFilesPath + " ./histo.tmp";
    Process cmdProc = Runtime.getRuntime().exec(command);
    cmdProc.waitFor();
    
	//System.out.println("main - numReducers : " + numReducers);
	long boundPortion = totalNbElements/numReducers+1;
	//System.out.println("main - totalNbElements : " + totalNbElements);
	//System.out.println("main - boundPortion : " + boundPortion);
	long currentCount = 0;
	int reducerID = 0;
	
	//System.out.println("main - reading histogram");
	
	//copy temp file in hdfs
	String histoFilePath = inputPath +"histo.tmp";
	command = "hdfs dfs -copyFromLocal ./histo.tmp "+ histoFilePath;
	cmdProc = Runtime.getRuntime().exec(command);
	cmdProc.waitFor();
	//System.out.println("exit value : " + cmdProc.exitValue());
	
	//fill treeMap - sorting anchors
	try{
        Path pt=new Path(histoFilePath);
        FileSystem fs = FileSystem.get(new Configuration());
        BufferedReader br=new BufferedReader(new InputStreamReader(fs.open(pt)));
        String line;
        line=br.readLine();
        while (line != null){
                //System.out.println(line);
                
                String tokens[] = line.split("\t");
                
                String anchor = tokens[0];
                long occurence = Long.parseLong(tokens[1]);
                
                histoTreeMap.put(anchor, occurence);
                //System.out.println("test - histoTreeMap.get(" + anchor + ") : " + histoTreeMap.get(anchor));
                line=br.readLine();
                
        }
	}catch(Exception e){
		System.out.println("\t\tERROR WHILE READING");
	}
	
	//System.out.println("\n\nmain - determine reducers");
	
	//determine reducersIDs
	
	switch (partitionVersion) {
	
		case 1: fillReducerTreeMapV1(totalNbElements, numReduceTasks);
		break;
		case 2: fillReducerTreeMapV2(totalNbElements, numReduceTasks);
		break;
		case 3: fillReducerTreeMapV3(totalNbElements, numReduceTasks);
		break;
		case 4: fillInfosReducersTreeMapV4(totalNbElements, numReduceTasks);
		break;
	}
	
	//OLD VERSION ALL ALGOS
	/*	
	
	long nbRemainingElements = totalNbElements;
	int nbRemainingReducers = numReducers;
	long limit = boundPortion;
	long distance = boundPortion;
	long newDistance = boundPortion;
	
	 for(Map.Entry<String,Long> entry : histoTreeMap.entrySet()) {
		String anchor = entry.getKey().toString();
		long occurence = entry.getValue();
		
		//System.out.println(anchor);
		//System.out.println(occurence);
         
        currentCount += occurence;
        //System.out.println("test - occurence : " + occurence);
        //System.out.println("test - currentCount : " + currentCount);
        //newDistance = Math.abs(boundPortion - currentCount);
        
        //System.out.println("test - newDistance : " + newDistance);
			  
			// version 3    
			//    if (newDistance <= distance)
			//   {
			// version 3 end  
			    	
				  reducerTreeMap.put(anchor, reducerID);
				  //distance = newDistance;
			  
			//version 3
			    //}
			    //else{
			    //	if (reducerID < numReduceTasks-1){
			    //		++reducerID;
			    //	}
			    //	reducerTreeMap.put(anchor, reducerID);
			    	
			    //	nbRemainingElements -= currentCount;
				//	nbRemainingReducers -= 1;
				//	if (nbRemainingReducers > 0){
				//		boundPortion = nbRemainingElements/nbRemainingReducers+1;
				//		//distance = boundPortion;
				//		System.out.println("test - new boundPortion : " + boundPortion);
				//	}
					
				//	currentCount = 0;
				//	distance = boundPortion;
				//	newDistance = 0;
			    //}
			// version 3 end

			//version 1,2	  
			if (currentCount > limit){

				++reducerID;
				  
				//reducerTreeMap.put(anchor, reducerID);
				
				//readjust bound portion
				
				//
				//version 1 
				//try to balance charge on remaining reducers
				//(the first one will probably be the more charged reducer)
				
				nbRemainingElements -= currentCount;
				nbRemainingReducers -= 1;
				boundPortion = nbRemainingElements/nbRemainingReducers+1;
				currentCount = 0;
				limit = boundPortion;
				
				
				//
				// version 2
				// balance by portions :
				// the first one will at least work on the first n ones
				// the second from n to 2n... etc
				// if one portion reaches the limit of the next reducer,
				// then the next one won't work at all...
				//
				//limit += boundPortion;
				//System.out.println("test - limit : " + limit);
				//
				
				//
				// version 3
				// take the next set in the reducer if the distance to
				// bound portion get smaller, and then reajust the perfect
				// bound
				// else stop
				//
				
			}
			//version 1,2 end
         
         //System.out.println("test - reducerTreeMap.get(" + anchor + ") : " + reducerTreeMap.get(anchor));
         //enter anchor and occurences in treeMap
         //keyReducer.put(key, reducerID);
         //count total nb occurences
         //totalNbElements += val;
	
         
 }*/
	
	//save the treeMap in a file
	String treeMapFilePath = inputPath +"treeMap.tmp";
	Path treeMapFile = new Path(treeMapFilePath);
	if ( hdfs.exists( treeMapFile )) { hdfs.delete( treeMapFile, true ); } 
	OutputStream os = hdfs.create( treeMapFile/*,
	    new Progressable() {
	        public void progress() {
	            out.println("...bytes written: [ "+bytesWritten+" ]");
	        } }*/);
	BufferedWriter bw = new BufferedWriter( new OutputStreamWriter( os, "UTF-8" ) );
	 
	if (partitionVersion < 4) {
		 
		 for(Map.Entry<String, Integer> entry : reducerTreeMap.entrySet()) {
			 
				String anchor = entry.getKey();
				int ID = entry.getValue();
	
				bw.write(anchor + ":" + ID + "\n");
			 }	
	}
	else {
		 
		for(Map.Entry<String, InfosReducers> entry : infosReducersTreeMap.entrySet()) {
			 
			String anchor = entry.getKey();
			InfosReducers ir = entry.getValue();
			int numFirstReducer = ir.getNumActualReducerID();
			long nbTasksLeftForActualReducer = ir.getNbTasksLeftForActualReducer();
	
			bw.write(anchor + ":" + numFirstReducer + ":" + nbTasksLeftForActualReducer +"\n");
		 }
	}
	
	 bw.close();
	 
	 /*String sreducerTreeMap = "";
	 
	 for(Map.Entry<String, Integer> entry : reducerTreeMap.entrySet()) {
		 
		String anchor = entry.getKey();
		int ID = entry.getValue();
		
		sreducerTreeMap += anchor + ":" + ID + ";";
	 }
	 sreducerTreeMap = sreducerTreeMap.substring(0, sreducerTreeMap.length() - 1);
	 System.out.println(sreducerTreeMap);*/
	 
	 conf.set("treeMapFilePath", treeMapFilePath);
	 conf.set("partitionVersion", Integer.toString(partitionVersion));
	 long limitReducerCapacity = totalNbElements/numReducers+1;
	 conf.set("limitReducerCapacity", Long.toString(limitReducerCapacity));
	
	//delete temp files
	command = "hdfs dfs -rm " + histoFilePath;
	cmdProc = Runtime.getRuntime().exec(command);
	cmdProc.waitFor();
	command = "rm ./histo.tmp";
	cmdProc = Runtime.getRuntime().exec(command);
	cmdProc.waitFor();
	//command = "hdfs dfs -rm " + histosFilesPath;
	command = "hdfs dfs -rm -r " + histosFilesPath;//hdfs://bilbo:9000/user/tbraquel/res/res2red.tmp";
	cmdProc = Runtime.getRuntime().exec(command);
	cmdProc.waitFor();
 
	//System.out.println("\n\nmain - sort map");
    
     //SORTMAP
     
	//numReduceTasks = 4;
    //Create Job
    Job job2 = Job.getInstance(conf, "Sort");  
    job2.setJarByClass(Sort.class);
    
    //Number of Reducer tasks.
    job2.setNumReduceTasks(numReduceTasks);
    numReducers = job2.getNumReduceTasks();
    //System.out.println(numReducers);
    
    // File Input and Output paths
    FileInputFormat.addInputPath(job2, new Path(args[0]));
    FileOutputFormat.setOutputPath(job2, new Path(args[1]));
    
    //Set Mapper class and Output format for key-value pair
    job2.setMapperClass(SortMapper.class);
    job2.setMapOutputKeyClass(Text.class);
    job2.setMapOutputValueClass(Text.class);
    
    //set partitioner statement
    job2.setPartitionerClass(SortPartitioner.class);
    
    //Set Reducer class and Input/Output format for key-value pair
    job2.setReducerClass(SortReducer.class);
    
    //Input and Output format for data
    job2.setInputFormatClass(TextInputFormat.class);
    job2.setOutputFormatClass(TextOutputFormat.class);
    job2.setOutputKeyClass(Text.class);
    job2.setOutputValueClass(Text.class);
    
    
    int exitMsg = job2.waitForCompletion(true) ? 0 : 1;
    
    //delete tmp files
    hdfs.delete(treeMapFile, false);
    hdfs.close();
    
    System.exit(exitMsg);
    
    
    
    /*
    // Create job and parse CLI parameters
    Job job = Job.getInstance(conf, "Sort");
    job.setJarByClass(Sort.class);
    
    Path inputPath = new Path(args[0]);
    Path partitionOutputPath = new Path(args[1]);
    Path outputPath = new Path(args[2]);

    // The following instructions should be executed before writing the partition file
    int numReduceTasks = 2;
    job.setNumReduceTasks(numReduceTasks);
    FileInputFormat.setInputPaths(job, inputPath);
    TotalOrderPartitioner.setPartitionFile(job.getConfiguration(), partitionOutputPath);
    job.setInputFormatClass(TextInputFormat.class);
    job.setMapOutputKeyClass(Text.class);
    job.setMapOutputValueClass(Text.class);
    job.setOutputFormatClass(TextOutputFormat.class);
    job.setOutputKeyClass(Text.class);
    job.setOutputValueClass(Text.class);
    
    // Write partition file with random sampler
    InputSampler.Sampler<Text,Text> sampler = new InputSampler.RandomSampler<>(0.01, numReduceTasks, numReduceTasks-1);
    InputSampler.writePartitionFile(job, sampler);

    // Use TotalOrderPartitioner and default identity mapper and reducer 
    job.setPartitionerClass(TotalOrderPartitioner.class);
    job.setMapperClass(SortMapper.class);
    job.setReducerClass(SortReducer.class);

    FileOutputFormat.setOutputPath(job, outputPath);
    System.exit(job.waitForCompletion(true) ? 0 : 1);
    */
    
  }
}
