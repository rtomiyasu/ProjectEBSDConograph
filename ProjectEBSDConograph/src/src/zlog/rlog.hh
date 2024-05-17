/*
 * The MIT License

   Conograph (powder auto-indexing program)

Copyright (c) <2012> <Ryoko Oishi-Tomiyasu, KEK>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
 *
 */
#ifndef RLOG_HH_
#define RLOG_HH_

#ifdef _OPENMP
	#include <omp.h>
#endif
#include <vector>
#include <string>

using namespace std;


typedef int zListnerID;

/**********************************************************************************/
/**
 * log class
 */
/**
 * level of error
 */
enum zLevel
{
    zDebug,                 ///< debug
    zInfo,                  ///< infomation
    zWarn,                  ///< Warning
    zError,                 ///< Error
    zFatal,                 ///< Fatal error
};


class CRLog
{
public:
    /** Log listener class to receive messages from everywhere in the program.
     *  Describe how to record log messages by inheritting this class.
     */
    class Listner
    {
    public:
       	/// Destructor.
       	virtual ~Listner() {};
       	/// Receive a message.
       	virtual void log(const zLevel& cLevel, ///< [in]:level of error
                         const string& cpStr ///< [in]: message of string type
                       ) = 0;

       	/// Receive a message.(The ID of the thread which calls this function is available here.)
       	virtual void log(const zLevel& cLevel, ///< [in]:level of error
       					 const int& threadID, ///< [in]:ID of thread
                         const string& cpStr ///< [in]:message of string type
                       ) = 0;

       	/** Returns class ID.
         * This ID is necessary remove this class from a list of registered listener.
         */
       	virtual zListnerID putId(void) const = 0;
    };

    /**
     * helper class to set a thread number and an error level.
     */
    class set
    {
    private:
        zLevel m_Level;
        int threadID;
    public:
        set(const zLevel& cLevel) : m_Level(cLevel), threadID(-1) {};
        set(const zLevel& cLevel, const int& tID) : m_Level(cLevel), threadID(tID) {};
        inline void log(const string& cpStr);
    };

    static bool append(Listner* cpPtr);
    static bool remove(const zListnerID& cId);
    static void clear(void);

    ~CRLog(void);
    static CRLog* getInstance(void);
    void log(const zLevel& cLevel, const string& cpStr);
    void log(const zLevel& cLevel, const int& threadID, const string& cpStr);

private:
    vector<Listner*> m_ListnerAry;

    CRLog(void);
    CRLog(const CRLog&);
    bool append_(Listner* cpPtr);
    bool remove_(const zListnerID& cId);
    void clear_(void);
};

/**
 * Send a message of string type.
 */
inline void CRLog::set::log(const string& cpStr)         ///< [in]:string to output
{
	if( threadID < 0 )
	{
		CRLog::getInstance()->log(m_Level, cpStr);
	}
	else
	{
		CRLog::getInstance()->log(m_Level, threadID, cpStr);
	}
};

#define ZLOG_INFO CRLog::set(zInfo).log
#define ZLOG_WARN CRLog::set(zWarn).log
#define ZLOG_ERROR CRLog::set(zError).log
#define ZLOG_FATAL CRLog::set(zFatal).log

#ifdef _OPENMP
#define ZLOG_INFO_THREAD CRLog::set(zInfo, omp_get_thread_num()).log
#define ZLOG_WARN_THREAD CRLog::set(zWarn, omp_get_thread_num()).log
#define ZLOG_ERROR_THREAD CRLog::set(zError, omp_get_thread_num()).log
#define ZLOG_FATAL_THREAD CRLog::set(zFatal, omp_get_thread_num()).log
#endif

#endif
