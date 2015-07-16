package database;

import java.util.Iterator;

import main.Thunder;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

public class ListDBTables {

	/**
	 * Specify the available command-line interface options
	 * @return
	 */
	@SuppressWarnings("static-access")
	public static Options getCmdLineOptions(){
		Options options = new Options();
		options.addOption(OptionBuilder.withArgName("databasePath").hasArg().withDescription("query the database at this path\n(e.g. /path/to/db)").create("d"));
		return options;
	}

	
	public static void main(String[] args) throws Exception {
		
		CommandLine cmdArgs = Thunder.parseArgs(args, getCmdLineOptions());
		if(cmdArgs.hasOption("d")){
			DBConnect db = new DBConnect_SQLite(cmdArgs.getOptionValue("d"));	
			Iterator<String> it = db.getTableNames().iterator();
			while(it.hasNext())
				System.out.println(it.next());
			db.closeConnection();
		}else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(Thunder.THUNDER_EXE_COMMAND+" ListDBTables", getCmdLineOptions());
			System.out.println();
		}
		
		
	}

}
