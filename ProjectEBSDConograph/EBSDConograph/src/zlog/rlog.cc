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
#include "../zlog/rlog.hh"

#include <assert.h>

/**
 * constructor
 */
CRLog::CRLog(void)
{
}

/**
 * copy constructor
 */
CRLog::CRLog(const CRLog&)
{
}

/**
 * destructor
 */
CRLog::~CRLog(void)
{
    clear();
}

/**
 * Gets the instance of singleton class.
 */
CRLog* CRLog::getInstance(void)
{
	static CRLog aPtr;
    return &aPtr;
};

/**
 * Registration of listner classes
 * @retval true -> success
 * @retval false -> fail
 * @attention If this method returns true, the pointer is deleted in this class.
 * If it returns false, the pointer is not deleted.
 */
bool CRLog::append(Listner* cpPtr)            ///< [in]:listener class to register.
{
	bool ans;
#ifdef _OPENMP
	#pragma omp critical
#endif
	{
		ans = getInstance()->append_(cpPtr);
	}
	return ans;
};

/**
 * Remove a listener class.
 * @retval true -> success
 * @retval false -> fail
 */
bool CRLog::remove(const zListnerID& cId)                    ///< [in]:ID of a listner class
{
	bool ans;
#ifdef _OPENMP
	#pragma omp critical
#endif
	{
		ans = getInstance()->remove_(cId);
	}
	return ans;
};

/**
 * Erase all the registered listner classes.
 */
void CRLog::clear(void)
{
#ifdef _OPENMP
	#pragma omp critical
#endif
	{
		getInstance()->clear_();
	}
};

/**
 * Send a message to every registered listener class
 */
void CRLog::log(const zLevel& cLevel,                 ///< [in]:level of error
    const string& cpStr)         ///< [in]:message of string type
{
#ifdef _OPENMP
	#pragma omp critical
#endif
	{
	    for(vector<Listner*>::iterator i=m_ListnerAry.begin(); i != m_ListnerAry.end(); ++i)
	    {
	        (*i)->log(cLevel, cpStr);
	    }
	}
};

/**
 * Send a message to every registered listener class
 */
void CRLog::log(const zLevel& cLevel,                 ///< [in]:level of error,
    const int& threadID,        ///< [in]:ID of the thread which calls this method,
    const string& cpStr)         ///< [in]:message of string type.
{
#ifdef _OPENMP
	#pragma omp critical
#endif
	{
	    for(vector<Listner*>::iterator i=m_ListnerAry.begin(); i != m_ListnerAry.end(); ++i)
	    {
	        (*i)->log(cLevel, threadID, cpStr);
	    }
	}
};

/**
 * Registration of Listener classes
 * @retval true -> success
 * @retval false -> fail
 * @see CRLog::append
 */
bool CRLog::append_(Listner* cpPtr)            ///< [in]:Instance to register
{
	const zListnerID aId = cpPtr->putId();
    for(vector<Listner*>::iterator i =m_ListnerAry.begin(); i != m_ListnerAry.end(); ++i)
    {
        if( (*i)->putId() == aId ) return false;
    }
    m_ListnerAry.push_back(cpPtr);
    return true;
};

/**
 * Remove of Listener classes
 * @retval true -> success
 * @retval false -> fail
 * @see CRLog::remove
 */
bool CRLog::remove_(const zListnerID& cId) ///< [in]:ID of a listener class
{
    for(vector<Listner*>::iterator i=m_ListnerAry.begin(); i!=m_ListnerAry.end(); ++i)
    {
        if( (*i)->putId() == cId )
        {
            delete *i;
            m_ListnerAry.erase(i);
            return true;
        }
    }
    return false;
};

/**
 * Erase of all the listener classes.
 * @retval true -> success
 * @retval false ->fail
 * @see CRLog::clear
 */
void CRLog::clear_(void)
{
    for(vector<Listner*>::iterator i=m_ListnerAry.begin();i!=m_ListnerAry.end(); ++i)
    {
        delete *i;
    }
    m_ListnerAry.clear();
};
