package database;

import java.io.File;
import java.sql.DriverManager;

public class DBConnect_SQLite extends DBConnect {
	
	/**
	 * 
	 * @param dbPath
	 * @throws Exception
	 */
	public DBConnect_SQLite(String dbPath) throws Exception {
		connect(dbPath, "", "");
	}
	
	
	/**
	 * 
	 * @param dbPath
	 * @param userName
	 * @param password
	 * @throws Exception
	 */
	public DBConnect_SQLite(String dbPath, String userName, String password) throws Exception {
		connect(dbPath, userName, password);
	}
	
	
	/**
	 * Obtain a connection to the database
	 * @param dbPath
	 * @param userName
	 * @param password
	 * @throws Exception
	 */
	protected void connect(String dbPath, String userName, String password) throws Exception {
		Class.forName("org.sqlite.JDBC");
		System.out.println("Connecting to DB at: '"+(new File(dbPath)).getAbsolutePath()+"'");
		this.setConnection(DriverManager.getConnection("jdbc:sqlite:"+(new File(dbPath)).getAbsolutePath(), userName, password));
		
		this.setQuery_ContainsTable("SELECT name FROM sqlite_master WHERE type='table'");
	}
	
	

	public static void main(String[] a) throws Exception {
		DBConnect_SQLite db = new DBConnect_SQLite("/Users/robk/Desktop/TESTDB.sqlite.db");
//		db.getStatement().executeQuery("SELECT...");
		db.closeConnection();
	}

}
