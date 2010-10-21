package nqcmol.tools;

import java.io.IOException;
import java.io.Writer;
import java.util.Stack;
import java.util.logging.Level;
import java.util.logging.Logger;


/**
 * Makes writing XML much much easier. 
 *
 * @author <a href="mailto:bayard@generationjava.com">Henri Yandell</a>
 * @version 0.1
 */
public class XmlWriter {

    private Writer writer;      // underlying writer
    private Stack stack;        // of xml entity names
    private StringBuffer attrs; // current attribute string
    private boolean empty;      // is the current node empty
    private boolean closed;     // is the current node closed...

    /**
     * Create an XmlWriter on top of an existing java.io.Writer.
     */
    public XmlWriter(Writer writer) {
        this.writer = writer;
        this.closed = true;
        this.stack = new Stack();
    }

    /**
     * Begin to output an entity. 
     *
     * @param String name of entity.
     */
    public XmlWriter writeEntity(String name)  {
        try {
            closeOpeningTag();
            this.closed = false;
            this.writer.write("<");
            this.writer.write(name);
            stack.add(name);
            this.empty = true;
            
        } catch (IOException ex) {
            Logger.getLogger(XmlWriter.class.getName()).log(Level.SEVERE, null, ex);
        }
        return this;
    }

    // close off the opening tag
    private void closeOpeningTag() {
        if (!this.closed) {
            try {
                writeAttributes();
                this.closed = true;
                this.writer.write(">\n");
            } catch (IOException ex) {
                Logger.getLogger(XmlWriter.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    // write out all current attributes
    private void writeAttributes() {
        if (this.attrs != null) {
            try {
                this.writer.write(this.attrs.toString());
                this.attrs.setLength(0);
                this.empty = false;
            } catch (IOException ex) {
                Logger.getLogger(XmlWriter.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    /**
     * Write an attribute out for the current entity. 
     * Any xml characters in the value are escaped.
     * Currently it does not actually throw the exception, but 
     * the api is set that way for future changes.
     *
     * @param String name of attribute.
     * @param String value of attribute.
     */
    public XmlWriter writeAttribute(String attr, String value)  {
        if (this.attrs == null) {
            this.attrs = new StringBuffer();
        }
        this.attrs.append(" ");
        this.attrs.append(attr);
        this.attrs.append("=\"");
        this.attrs.append(escapeXml(value));
        this.attrs.append("\"");
        return this;
    }

    /**
     * End the current entity. This will throw an exception 
     * if it is called when there is not a currently open 
     * entity.
     */
    public XmlWriter endEntity(){
        try {
            if(this.stack.empty()) {

                throw new IOException("Called endEntity too many times. ");

            }
            String name = (String)this.stack.pop();
            if (name != null) {
                if (this.empty) {
                    writeAttributes();
                    this.writer.write("/>");
					this.writer.write("\n");
                } else {
                    this.writer.write("</");
                    this.writer.write(name);
                    this.writer.write(">");
					this.writer.write("\n");
                }
                this.empty = false;
                this.closed = true;
            }
        } catch (IOException ex) {
            Logger.getLogger(XmlWriter.class.getName()).log(Level.SEVERE, null, ex);
        }
        return this;
    }

    /**
     * Close this writer. It does not close the underlying 
     * writer, but does throw an exception if there are 
     * as yet unclosed tags.
     */
    public void close() throws IOException {
        if(!this.stack.empty()) {
            throw new IOException("Tags are not all closed. "+
                "Possibly, "+this.stack.pop()+" is unclosed. ");
        }
    }

    /**
     * Output body text. Any xml characters are escaped. 
     */
    public XmlWriter writeText(String text)  {
        try {
            closeOpeningTag();
            this.empty = false;
            this.writer.write(escapeXml(text));            
        } catch (IOException ex) {
            Logger.getLogger(XmlWriter.class.getName()).log(Level.SEVERE, null, ex);
        }
        return this;
    }

	/**
     * Output body text
     */
    public XmlWriter writeNormalText(String text)  {
        try {
            closeOpeningTag();
            this.empty = false;
            this.writer.write(text);
           
        } catch (IOException ex) {
            Logger.getLogger(XmlWriter.class.getName()).log(Level.SEVERE, null, ex);
        }
         return this;
    }

	 public XmlWriter flush(){
        try {
            this.writer.flush();
        } catch (IOException ex) {
            Logger.getLogger(XmlWriter.class.getName()).log(Level.SEVERE, null, ex);
        }
        return this;
    }

    // Static functions lifted from generationjava helper classes
    // to make the jar smaller.
    
    // from XmlW
   
    static public String escapeXml(String str) {
        str = replaceString(str,"&","&amp;");
        str = replaceString(str,"<","&lt;");
        str = replaceString(str,">","&gt;");
        str = replaceString(str,"\"","&quot;");
        str = replaceString(str,"'","&apos;");
        return str;
    }  
    

    // from StringW
    /*
    static public String replaceString(String text, String repl, String with) {
        return replaceString(text, repl, with, -1);
    }  
    */
    /**
     * Replace a string with another string inside a larger string, for
     * the first n values of the search string.
     *
     * @param text String to do search and replace in
     * @param repl String to search for
     * @param with String to replace with
     * @param n    int    values to replace
     *
     * @return String with n values replacEd
     */
     
    static public String replaceString(String text, String repl, String with) {
        if(text == null) {
            return null;
        }
 
        StringBuffer buffer = new StringBuffer(text.length());
        int start = 0;
        int end = 0;
        while( (end = text.indexOf(repl, start)) != -1 ) {
            buffer.append(text.substring(start, end)).append(with);
            start = end + repl.length();
 
        }
        buffer.append(text.substring(start));
 
        return buffer.toString();
    }              
    

    static public void test1() throws IOException {
        Writer writer = new java.io.StringWriter();
        XmlWriter xmlwriter = new XmlWriter(writer);
        xmlwriter.writeEntity("person").writeAttribute("name", "fred").writeAttribute("age", "12").writeEntity("phone").writeText("4254343").endEntity().writeEntity("bob").endEntity().endEntity();
        xmlwriter.close();
        System.err.println(writer.toString());
    }
    static public void test2() throws IOException {
        Writer writer = new java.io.StringWriter();
        XmlWriter xmlwriter = new XmlWriter(writer);
        xmlwriter.writeEntity("person");
        xmlwriter.writeAttribute("name", "fred");
        xmlwriter.writeAttribute("age", "12");
        xmlwriter.writeEntity("phone");
        xmlwriter.writeText("4254343");
        xmlwriter.endEntity();
        xmlwriter.writeEntity("bob");
        xmlwriter.endEntity();
        xmlwriter.endEntity();
        xmlwriter.close();
        System.err.println(writer.toString());
    }

}
