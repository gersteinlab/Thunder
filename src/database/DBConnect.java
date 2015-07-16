package database;

import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;

public abstract class DBConnect {
	
	private Connection connection;
	private String query_containsTable = "";
	
	protected abstract void connect(String dbPath, String userName, String password) throws Exception;
	
	
	
	/**
	 * Setters
	 * @throws SQLException 
	 */
	protected void setConnection(Connection conn) throws SQLException{ this.connection = conn; }
	protected void setQuery_ContainsTable(String sql){ this.query_containsTable = sql; }
	
	/**
	 * Getters
	 */
//	private Connection getConnection(){ return this.connection; }
//	protected void getQuery_ContainsTable(String sql){ this.query_containsTable = sql; }
	
	/**
	 * Return the connection to the database
	 * @throws SQLException 
	 */
	public void closeConnection() throws SQLException{ this.connection.close(); }
	
	/**
	 * Set whether to auto commit changes to the DB
	 * @param autoCommit
	 * @throws SQLException
	 */
	public void setAutoCommit(boolean autoCommit) throws SQLException{ this.connection.setAutoCommit(autoCommit); }
	
	/**
	 * Commit the SQL command stack to the database
	 * @throws SQLException
	 */
	public void commit() throws SQLException { this.connection.commit(); }
	
	
	/**
	 * Return a Statement object for this database
	 * @throws SQLException 
	 */
	public Statement createStatement() throws SQLException{ return this.connection.createStatement(); }
	
	
	/**
	 * Return a PreparedStatement object for this database
	 * @throws SQLException 
	 */
	public PreparedStatement getPreparedStatement(String sql) throws SQLException{
		return this.connection.prepareStatement(sql);
	}
	
	
	/**
	 * Check this DB to see if this table exists
	 */
	public boolean containsTable(String tableName) throws SQLException {
		return getTableNames().contains(tableName);
	}
	

	/**
	 * 
	 * @return
	 * @throws SQLException
	 */
	public ArrayList<String> getTableNames() throws SQLException {
		ArrayList<String> names = new ArrayList<String>();
		
		Statement st = this.connection.createStatement();
		ResultSet rs = st.executeQuery(this.query_containsTable);	
		while(rs.next()){
			names.add(rs.getString(1));
		}
		st.close();
		return names;
	}
	
}
