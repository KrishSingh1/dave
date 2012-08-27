///////////////////////////////////////////////////////////////////////////////////////
/// \file archive.h
/// \brief Classes to make (de)serializing to/from streams convenient
///
/// $Date$
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_ARCHIVE_H
#define LPJ_GUESS_ARCHIVE_H

#include <ostream>
#include <istream>

/// Abstract base class for ArchiveInStream and ArchiveOutStream
/** The base class declares the transfer function, which will read
 *  data from a stream in ArchiveInStream, and write data to a stream
 *  in ArchiveOutStream. By having the same interface for both
 *  cases, classes can have one function for both serializing and
 *  deserializing.
 *
 *  Sometimes we do need to know which direction the ArchiveStream
 *  is working in though, so that can be queried with the save
 *  function.
 */
class ArchiveStream {
public:

	/// Checks if this ArchiveStream is saving data to stream or reading
	virtual bool save() const = 0;

	/// Write or read data to/from a stream
	/** \param s   Data buffer to write out, or read to.
	 *             Needs to be n bytes.
	 *  \param n   Number of bytes to read or write
	 */
	virtual void transfer(char* s, std::streamsize n) = 0;
};

/// Class for reading data from an istream
/** \see ArchiveStream for more documentation */
class ArchiveInStream : public ArchiveStream {
public:
	ArchiveInStream(std::istream& strm);

	bool save() const;

	void transfer(char* s, std::streamsize n);

private:
	/// The stream we're reading from
	std::istream& in;
};

/// Class for writing data to an ostream
/** \see ArchiveStream for more documentation */
class ArchiveOutStream : public ArchiveStream {
public:
	ArchiveOutStream(std::ostream& strm);

	bool save() const;

	void transfer(char* s, std::streamsize n);

private:
	/// The stream we're writing to
	std::ostream& out;
};

/// Interface showing that a class can serialize itself
/** Classes which can serialize themselves through an ArchiveStream
 *  should inherit from this class, and implement the serialize function.
 */
class Serializable {
public:
	/// Needs to be implemented by all sub-classes
	virtual void serialize(ArchiveStream& arch) = 0;
};

/// Operator overloading of &, allowing serialization to be chained
/** Since the operator returns the original stream, we can serialize
 *  multiple items in the following way:
 *
 *  \code
 *    ArchiveOutStream arch(os);
 *
 *    arch & width & height & name;
 *  \endcode
 *
 *  This function can be used for any item where the memory representation
 *  of the item can be written to file bit by bit, such as the primitive
 *  data types. Classes should typically do their own serialization by
 *  implementing Serializable.
 */
template<typename T>
ArchiveStream& operator&(ArchiveStream& stream, T& data) {
	stream.transfer((char*)&data, sizeof(data));
	return stream;
}

/// Operator overloading of &, for classes implementing Serializable
/** This makes it possible to serialize complex objects with the &-syntax,
 *  given that they implement Serializable.
 */
ArchiveStream& operator&(ArchiveStream& stream, Serializable& obj);

#endif // LPJ_GUESS_ARCHIVE_H
