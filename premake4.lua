#!lua

solution "stlcut"
	location ( "." )
	targetdir("build")
	configurations { "Debug", "Release" }
	platforms {"native", "x64", "x32"}

	configuration "Debug"
		defines { "DEBUG" }
		flags { "Symbols", "ExtraWarnings"}

	configuration "Release"
		defines { "NDEBUG" }
		flags { "Optimize", "ExtraWarnings"}    


	project "lib"
		language "C++"
		kind "SharedLib"
		files { "stlcut.cpp", "*.h" }
		targetname ("stlcut")
		links { "admesh", "poly2tri" }
		configuration { "linux" }
			targetextension (".so.1")
			linkoptions { "-Wl,-soname,libstlcut.so.1" }
			postbuildcommands { "ln -sf libstlcut.so.1 build/libstlcut.so" }
			--cleancommands { "rm build/libstlcut.so || :" }
		configuration { "macosx" }
			targetextension (".1.dylib")
			postbuildcommands { "ln -sf libstlcut.1.dylib build/libstlcut.dylib" }
			--cleancommands { "rm build/libstlcut.dylib || :" }

	project "cut"
		kind "ConsoleApp"
		language "C++"
		linkoptions { "-Lbuild" }
		files { "prgstlcut.cpp" }
		targetname ("stlcut")
		links { "lib", "admesh", "poly2tri" }
